#!/usr/bin/env python3
"""Flask web front-end for multi-genome ANI/AAI comparison.

Run:
    .venv/bin/python app.py
Then open http://localhost:5000 in your browser.
"""

from __future__ import annotations

import base64
import json
import math
import os
import shutil
import sqlite3
import subprocess
import tempfile
import uuid
from datetime import datetime, timezone
from pathlib import Path

from flask import Flask, jsonify, render_template, request

import aai_calculation
import generate_heatmap

app = Flask(__name__)
app.config["MAX_CONTENT_LENGTH"] = 500 * 1024 * 1024  # 500 MB upload limit

PROJECT_DIR = Path(__file__).parent.resolve()
PYANI_BIN = PROJECT_DIR / ".venv" / "bin" / "pyani-plus"
OUTPUTS_DIR = PROJECT_DIR / "outputs"
FASTA_EXTENSIONS = {".fa", ".fas", ".fasta", ".fna"}

ANI_RESULTS_DB = PROJECT_DIR / "ani_results.db"
AAI_RESULTS_DB = PROJECT_DIR / "aai_results.db"
RUN_CACHE_LIMIT = 10

OUTPUTS_DIR.mkdir(exist_ok=True)


def _is_fasta(filename: str) -> bool:
    return Path(filename).suffix.lower() in FASTA_EXTENSIONS


def _utc_now_iso() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def _db_connection(db_path: Path) -> sqlite3.Connection:
    conn = sqlite3.connect(db_path)
    conn.row_factory = sqlite3.Row
    return conn


def _results_db_for(analysis_type: str) -> Path:
    if analysis_type == "ani":
        return ANI_RESULTS_DB
    if analysis_type == "aai":
        return AAI_RESULTS_DB
    raise ValueError(f"Unknown analysis type: {analysis_type}")


def _init_results_db() -> None:
    for db_path in (ANI_RESULTS_DB, AAI_RESULTS_DB):
        with _db_connection(db_path) as conn:
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS run_cache (
                    run_id TEXT PRIMARY KEY,
                    analysis_type TEXT NOT NULL CHECK (analysis_type IN ('ani', 'aai')),
                    created_at TEXT NOT NULL,
                    last_accessed_at TEXT NOT NULL,
                    image_path TEXT NOT NULL,
                    table_json TEXT NOT NULL
                )
                """
            )
            conn.execute(
                """
                CREATE INDEX IF NOT EXISTS idx_run_cache_access
                ON run_cache(last_accessed_at DESC, created_at DESC)
                """
            )


def _evict_lru_runs(conn: sqlite3.Connection) -> None:
    rows = conn.execute(
        """
        SELECT run_id, image_path
        FROM run_cache
        ORDER BY last_accessed_at DESC, created_at DESC
        """
    ).fetchall()

    for row in rows[RUN_CACHE_LIMIT:]:
        conn.execute("DELETE FROM run_cache WHERE run_id = ?", (row["run_id"],))
        img = Path(row["image_path"])
        if img.exists():
            try:
                img.unlink()
            except OSError:
                pass


def _cache_run(analysis_type: str, run_id: str, heatmap_path: Path, table: dict) -> None:
    db_path = _results_db_for(analysis_type)
    now = _utc_now_iso()
    with _db_connection(db_path) as conn:
        conn.execute(
            """
            INSERT OR REPLACE INTO run_cache
              (run_id, analysis_type, created_at, last_accessed_at, image_path, table_json)
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            (
                run_id,
                analysis_type,
                now,
                now,
                str(heatmap_path.resolve()),
                json.dumps(table, separators=(",", ":")),
            ),
        )
        _evict_lru_runs(conn)


def _touch_run(analysis_type: str, run_id: str) -> None:
    db_path = _results_db_for(analysis_type)
    with _db_connection(db_path) as conn:
        conn.execute(
            "UPDATE run_cache SET last_accessed_at = ? WHERE run_id = ?",
            (_utc_now_iso(), run_id),
        )


_init_results_db()


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/history/<analysis_type>", methods=["GET"])
def history_list(analysis_type: str):
    if analysis_type not in {"ani", "aai"}:
        return jsonify(error="Unknown analysis type."), 400

    db_path = _results_db_for(analysis_type)
    with _db_connection(db_path) as conn:
        rows = conn.execute(
            """
            SELECT run_id, analysis_type, created_at, last_accessed_at
            FROM run_cache
            ORDER BY last_accessed_at DESC, created_at DESC
            LIMIT ?
            """,
            (RUN_CACHE_LIMIT,),
        ).fetchall()

    return jsonify(
        runs=[
            {
                "run_id": row["run_id"],
                "analysis_type": row["analysis_type"],
                "created_at": row["created_at"],
                "last_accessed_at": row["last_accessed_at"],
            }
            for row in rows
        ]
    )


@app.route("/history/run/<run_id>", methods=["GET"])
def history_run(run_id: str):
    # Find run in ANI cache first, then AAI cache.
    for analysis_type in ("ani", "aai"):
        db_path = _results_db_for(analysis_type)
        with _db_connection(db_path) as conn:
            row = conn.execute(
                """
                SELECT run_id, analysis_type, image_path, table_json
                FROM run_cache
                WHERE run_id = ?
                """,
                (run_id,),
            ).fetchone()

        if row is None:
            continue

        image_path = Path(row["image_path"])
        if not image_path.exists():
            with _db_connection(db_path) as conn:
                conn.execute("DELETE FROM run_cache WHERE run_id = ?", (run_id,))
            return jsonify(error="Cached image no longer exists."), 404

        _touch_run(analysis_type, run_id)
        img_b64 = base64.b64encode(image_path.read_bytes()).decode()
        table = json.loads(row["table_json"])

        return jsonify(
            run_id=row["run_id"],
            analysis_type=row["analysis_type"],
            image=img_b64,
            table=table,
        )

    return jsonify(error="Run not found in cache."), 404


@app.route("/analyze", methods=["POST"])
def analyze():
    uploaded = request.files.getlist("genomes")
    valid_files = [f for f in uploaded if f and f.filename and _is_fasta(f.filename)]

    if len(valid_files) < 2:
        return jsonify(error="Please upload at least 2 FASTA files (.fa / .fas / .fasta / .fna)."), 400

    run_id = uuid.uuid4().hex[:8]
    work_dir = Path(tempfile.mkdtemp(prefix=f"ani_run_{run_id}_"))

    try:
        for f in valid_files:
            dest = work_dir / Path(f.filename).name
            f.save(dest)

        db_path = work_dir / "ani_results.db"
        heatmap_path = OUTPUTS_DIR / f"heatmap_{run_id}.png"

        conda_ani_bin = Path.home() / "miniconda3" / "envs" / "ani-tools" / "bin"
        env = os.environ.copy()
        path_parts = [str(work_dir), str(PYANI_BIN.parent)]
        if conda_ani_bin.exists():
            path_parts.append(str(conda_ani_bin))
        if env.get("PATH"):
            path_parts.append(env["PATH"])
        env["PATH"] = os.pathsep.join(path_parts)

        shim = work_dir / ".pyani-plus-private-cli"
        shim.write_text(
            "#!/usr/bin/env python3\n"
            "from pyani_plus.private_cli import app\n"
            "if __name__ == '__main__':\n"
            "    app()\n",
            encoding="utf-8",
        )
        shim.chmod(0o755)

        cmd = [
            str(PYANI_BIN),
            "fastani",
            str(work_dir),
            "--database",
            str(db_path),
            "--create-db",
        ]

        result = subprocess.run(cmd, capture_output=True, text=True, env=env)

        if result.returncode != 0:
            error_text = result.stderr or result.stdout or "Unknown error from pyani-plus."
            return jsonify(error=f"ANI calculation failed:\n{error_text}"), 500

        genome_names = [Path(f.filename).stem for f in valid_files]
        title = f"ANI Identity (%) - {', '.join(genome_names)}"
        generate_heatmap.generate(db_path, heatmap_path, title=title)

        img_b64 = base64.b64encode(heatmap_path.read_bytes()).decode()

        df = generate_heatmap.load_identity_matrix(db_path)
        table = {
            "columns": list(df.columns),
            "index": list(df.index),
            "data": [[None if (isinstance(v, float) and math.isnan(v)) else float(v) for v in row] for row in df.values.tolist()],
        }

        _cache_run("ani", run_id, heatmap_path, table)
        return jsonify(image=img_b64, table=table, run_id=run_id)

    except Exception as exc:
        return jsonify(error=str(exc)), 500

    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


@app.route("/analyze-aai", methods=["POST"])
def analyze_aai():
    """Accept protein FASTA uploads, run FastAAI all-vs-all, return heatmap + table."""
    uploaded = request.files.getlist("genomes")
    valid_files = [f for f in uploaded if f and f.filename and _is_fasta(f.filename)]

    if len(valid_files) < 2:
        return jsonify(error="Please upload at least 2 protein FASTA files (.faa / .fasta / .fa / .fas)."), 400

    run_id = uuid.uuid4().hex[:8]
    work_dir = Path(tempfile.mkdtemp(prefix=f"aai_run_{run_id}_"))

    try:
        staged: list[Path] = []
        for f in valid_files:
            dest = work_dir / Path(f.filename).name
            f.save(dest)
            staged.append(dest)

        heatmap_path = OUTPUTS_DIR / f"aai_heatmap_{run_id}.png"

        df = aai_calculation.compute_aai_matrix(staged)

        genome_names = [Path(f.filename).stem for f in valid_files]
        title = f"AAI Identity (%) - {', '.join(genome_names)}"
        generate_heatmap.generate_from_df(df, heatmap_path, title=title)

        img_b64 = base64.b64encode(heatmap_path.read_bytes()).decode()

        table = {
            "columns": list(df.columns),
            "index": list(df.index),
            "data": [[None if (isinstance(v, float) and math.isnan(v)) else float(v) for v in row] for row in df.values.tolist()],
        }

        _cache_run("aai", run_id, heatmap_path, table)
        return jsonify(image=img_b64, table=table, run_id=run_id)

    except Exception as exc:
        return jsonify(error=str(exc)), 500

    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


if __name__ == "__main__":
    port = int(os.environ.get("PORT", 8080))
    app.run(host="0.0.0.0", port=port, debug=True)
