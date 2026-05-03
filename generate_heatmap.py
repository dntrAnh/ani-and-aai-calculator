#!/usr/bin/env python3
"""Generate ANI heatmap images from a pyANI-plus SQLite database.

Usage (standalone):
    python generate_heatmap.py ani_results.db --run-id latest --output outputs/heatmap.png
"""

from __future__ import annotations

import argparse
import json
import sqlite3
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # non-interactive backend — must be set before pyplot import
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns


EXCEL_NUMBER_FORMAT = "0.############################"


def load_identity_matrix(db_path: Path, run_id: int | None = None) -> pd.DataFrame:
    """Return a labelled ANI identity DataFrame for the requested run.

    Args:
        db_path:  Path to the pyANI-plus SQLite database.
        run_id:   Specific run to load.  Defaults to the latest completed run.

    Returns:
        DataFrame with genome names as both index and columns,
        values as ANI percentage (0-100 scale).

    Raises:
        ValueError: If no completed run is found in the database.
    """
    conn = sqlite3.connect(db_path)
    try:
        cur = conn.cursor()

        if run_id is None:
            cur.execute(
                "SELECT run_id FROM runs WHERE status='Done' ORDER BY run_id DESC LIMIT 1"
            )
            row = cur.fetchone()
            if row is None:
                raise ValueError(
                    f"No completed ANI run found in {db_path}. "
                    "Run the analysis first."
                )
            run_id = row[0]

        # Load the identity matrix cached by pyANI-plus.
        cur.execute("SELECT df_identity FROM runs WHERE run_id=?", (run_id,))
        row = cur.fetchone()
        if row is None or row[0] is None:
            raise ValueError(f"Run {run_id} has no identity matrix stored.")

        df_json = json.loads(row[0])
        df = pd.DataFrame(
            data=df_json["data"],
            index=df_json["index"],
            columns=df_json["columns"],
        )

        # Map genome hashes -> friendly filenames.
        cur.execute(
            "SELECT genome_hash, fasta_filename FROM runs_genomes WHERE run_id=?",
            (run_id,),
        )
        hash_to_name = {
            h: Path(fn).stem for h, fn in cur.fetchall()
        }

        df.index = [hash_to_name.get(h, h) for h in df.index]
        df.columns = [hash_to_name.get(h, h) for h in df.columns]

        # Convert to percentage.
        if df.max().max() <= 1.0:
            df = df * 100

        return df
    finally:
        conn.close()


def plot_heatmap(
    df: pd.DataFrame,
    output_path: Path,
    title: str = "Identity (%)",
    vmin: float = 80,
    vmax: float = 100,
    cbar_label: str = "Identity (%)",
) -> Path:
    """Render a seaborn heatmap of an identity matrix and save it as PNG.

    Args:
        df:           Square identity DataFrame (values 0-100).
        output_path:  Destination PNG file.
        title:        Plot title.
        vmin:         Colour-scale minimum (default 80 for ANI, use 0 for AAI).
        vmax:         Colour-scale maximum (default 100).
        cbar_label:   Colour-bar label.

    Returns:
        The resolved path to the saved image.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    n = len(df)
    fig_size = max(5, n + 2)
    fig, ax = plt.subplots(figsize=(fig_size, fig_size))

    sns.heatmap(
        df,
        ax=ax,
        annot=True,
        fmt=".2f",
        cmap="YlOrRd",
        vmin=vmin,
        vmax=vmax,
        linewidths=0.5,
        linecolor="white",
        square=True,
        cbar_kws={"label": cbar_label, "shrink": 0.75},
    )

    ax.set_title(title, fontsize=14, fontweight="bold", pad=12)
    ax.set_xlabel("Reference Genome", fontsize=11)
    ax.set_ylabel("Query Genome", fontsize=11)
    ax.tick_params(axis="x", rotation=30, labelsize=9)
    ax.tick_params(axis="y", rotation=0, labelsize=9)

    plt.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

    return output_path.resolve()


def generate(
    db_path: Path,
    output_path: Path,
    run_id: int | None = None,
    title: str = "ANI Identity (%)",
) -> Path:
    """High-level helper: load ANI DB -> build heatmap -> save image.

    Returns the resolved path of the saved image.
    """
    df = load_identity_matrix(db_path, run_id=run_id)
    return plot_heatmap(df, output_path, title=title, vmin=80, vmax=100,
                        cbar_label="ANI (%)")


def generate_from_df(
    df: pd.DataFrame,
    output_path: Path,
    title: str = "AAI Identity (%)",
    vmin: float = 0,
    vmax: float = 100,
    cbar_label: str = "AAI (%)",
) -> Path:
    """Build a heatmap directly from a DataFrame (e.g. from FastAAI output).

    Returns the resolved path of the saved image.
    """
    return plot_heatmap(df, output_path, title=title, vmin=vmin, vmax=vmax,
                        cbar_label=cbar_label)


def export_df_to_excel(
    df: pd.DataFrame,
    output_path: Path,
    sheet_name: str = "Identity Matrix",
    index_label: str = "Query / Reference",
) -> Path:
    """Write an identity matrix to an Excel workbook without heatmap styling.

    Values are written as raw numeric cells and given a non-rounding display format
    so Excel shows the full stored value rather than forcing two decimals.
    """
    try:
        from openpyxl import load_workbook
    except ImportError as exc:
        raise RuntimeError(
            "Excel export requires the 'openpyxl' package. Install it with 'pip install openpyxl'."
        ) from exc

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with pd.ExcelWriter(output_path, engine="openpyxl") as writer:
        df.to_excel(writer, sheet_name=sheet_name, index_label=index_label)

    workbook = load_workbook(output_path)
    worksheet = workbook[sheet_name]
    worksheet.freeze_panes = "B2"

    for row in worksheet.iter_rows(min_row=2, min_col=2):
        for cell in row:
            if isinstance(cell.value, (int, float)):
                cell.number_format = EXCEL_NUMBER_FORMAT

    workbook.save(output_path)
    return output_path.resolve()


# ---------------------------------------------------------------------------
# CLI entry-point (standalone use)
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Generate ANI heatmap from pyANI-plus database.")
    p.add_argument("database", type=Path, help="Path to ani_results.db")
    p.add_argument(
        "-o", "--output",
        type=Path,
        default=Path("outputs/heatmap.png"),
        help="Output image path (default: outputs/heatmap.png)",
    )
    p.add_argument(
        "--run-id",
        type=int,
        default=None,
        help="Specific run ID to plot (default: latest completed run)",
    )
    p.add_argument(
        "--title",
        default="ANI Identity (%)",
        help="Heatmap title",
    )
    return p.parse_args()


if __name__ == "__main__":
    args = _parse_args()
    saved = generate(args.database, args.output, run_id=args.run_id, title=args.title)
    print(f"Heatmap saved: {saved}")
