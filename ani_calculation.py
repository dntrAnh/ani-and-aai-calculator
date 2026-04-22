#!/usr/bin/env python3
"""Basic wrapper to run ANI calculations with pyANI-plus.

Examples:
  python ani_calculation.py genomes_dir/
  python ani_calculation.py genomeA.fasta genomeB.fna -m fastani -d ani_results.db
"""

from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

import generate_heatmap


FASTA_EXTENSIONS = {".fa", ".fas", ".fasta", ".fna"}
DEFAULT_PYANI_BIN = str(Path(".venv/bin/pyani-plus")) if Path(".venv/bin/pyani-plus").exists() else "pyani-plus"


def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(
		description="Run ANI calculations with pyANI-plus from FASTA files/directories.",
	)
	parser.add_argument(
		"inputs",
		nargs="+",
		type=Path,
		help="One or more FASTA files or directories containing FASTA files.",
	)
	parser.add_argument(
		"-m",
		"--method",
		choices=["fastani", "anib", "anim", "dnadiff", "sourmash"],
		default="fastani",
		help="ANI method to run (default: fastani).",
	)
	parser.add_argument(
		"-d",
		"--database",
		type=Path,
		default=Path("ani_results.db"),
		help="SQLite database path for pyANI-plus runs (default: ani_results.db).",
	)
	parser.add_argument(
		"--run-name",
		default=None,
		help="Optional run name shown in pyANI-plus.",
	)
	parser.add_argument(
		"--executor",
		choices=["local", "slurm"],
		default="local",
		help="Execution backend for pyANI-plus (default: local).",
	)
	parser.add_argument(
		"--pyani-bin",
		default=DEFAULT_PYANI_BIN,
		help="Path to pyani-plus executable (default: auto-detect .venv/bin/pyani-plus, else pyani-plus).",
	)
	parser.add_argument(
		"--no-create-db",
		action="store_true",
		help="Do not pass --create-db, even if database does not exist.",
	)
	return parser.parse_args()


def is_fasta_file(path: Path) -> bool:
	return path.is_file() and path.suffix.lower() in FASTA_EXTENSIONS


def collect_fasta_files(inputs: list[Path]) -> list[Path]:
	fasta_files: list[Path] = []

	for input_path in inputs:
		if not input_path.exists():
			raise FileNotFoundError(f"Input does not exist: {input_path}")

		if input_path.is_dir():
			for child in sorted(input_path.iterdir()):
				if is_fasta_file(child):
					fasta_files.append(child.resolve())
		elif is_fasta_file(input_path):
			fasta_files.append(input_path.resolve())
		else:
			raise ValueError(
				f"Not a supported FASTA file or directory: {input_path} "
				f"(extensions: {', '.join(sorted(FASTA_EXTENSIONS))})"
			)

	# Keep first-seen order while deduplicating.
	unique_files = list(dict.fromkeys(fasta_files))
	return unique_files


def stage_fasta_files(fasta_files: list[Path], stage_dir: Path) -> None:
	used_names: set[str] = set()
	for src in fasta_files:
		target_name = src.name
		if target_name in used_names:
			stem, suffix = src.stem, src.suffix
			i = 1
			while f"{stem}_{i}{suffix}" in used_names:
				i += 1
			target_name = f"{stem}_{i}{suffix}"

		used_names.add(target_name)
		shutil.copy2(src, stage_dir / target_name)


def write_private_cli_shim(stage_dir: Path) -> Path:
	"""Create shim for pyANI-plus private CLI expected by some workflows."""
	shim = stage_dir / ".pyani-plus-private-cli"
	shim.write_text(
		"#!/usr/bin/env python3\n"
		"from pyani_plus.private_cli import app\n"
		"if __name__ == '__main__':\n"
		"    app()\n",
		encoding="utf-8",
	)
	shim.chmod(0o755)
	return shim


def run_pyani_plus(args: argparse.Namespace, staged_input_dir: Path) -> int:
	cmd = [
		args.pyani_bin,
		args.method,
		str(staged_input_dir),
		"--database",
		str(args.database),
		"--executor",
		args.executor,
	]

	if args.run_name:
		cmd.extend(["--name", args.run_name])

	if not args.no_create_db and not args.database.exists():
		cmd.append("--create-db")

	print("Running:")
	print(" ".join(cmd))
	print()

	# Ensure pyANI and ANI backend binaries are discoverable.
	env = os.environ.copy()
	path_parts = [str(staged_input_dir)]

	pyani_bin_path = Path(args.pyani_bin)
	if pyani_bin_path.parent != Path("."):
		path_parts.append(str(pyani_bin_path.expanduser().resolve().parent))

	conda_ani_bin = Path.home() / "miniconda3" / "envs" / "ani-tools" / "bin"
	if conda_ani_bin.exists():
		path_parts.append(str(conda_ani_bin))

	if env.get("PATH"):
		path_parts.append(env["PATH"])
	env["PATH"] = os.pathsep.join(path_parts)

	write_private_cli_shim(staged_input_dir)

	completed = subprocess.run(
		cmd,
		check=False,
		text=True,
		capture_output=True,
		env=env,
	)

	# Preserve pyANI-plus output for user visibility.
	if completed.stdout:
		print(completed.stdout, end="")
	if completed.stderr:
		print(completed.stderr, end="", file=sys.stderr)

	combined_output = f"{completed.stdout}\n{completed.stderr}"
	error_markers = ("Unhandled exception", "RuntimeError:", "Traceback")
	if completed.returncode == 0 and any(marker in combined_output for marker in error_markers):
		return 1

	return completed.returncode


def main() -> int:
	args = parse_args()

	try:
		fasta_files = collect_fasta_files(args.inputs)
	except (FileNotFoundError, ValueError) as exc:
		print(f"Error: {exc}", file=sys.stderr)
		return 2

	if len(fasta_files) < 2:
		print(
			"Error: ANI calculation needs at least 2 FASTA files. "
			"Provide multiple files or a directory with 2+ genomes.",
			file=sys.stderr,
		)
		return 2

	args.database = args.database.expanduser().resolve()
	args.database.parent.mkdir(parents=True, exist_ok=True)

	with tempfile.TemporaryDirectory(prefix="pyani_plus_input_") as tmp_dir:
		staged_input_dir = Path(tmp_dir)
		stage_fasta_files(fasta_files, staged_input_dir)

		return_code = run_pyani_plus(args, staged_input_dir)

	if return_code != 0:
		print("pyANI-plus run failed.", file=sys.stderr)
		return return_code

	print("pyANI-plus run completed successfully.")

	# Generate heatmap image — the sole user-facing output.
	output_dir = args.database.parent / "outputs"
	heatmap_path = output_dir / "ani_heatmap.png"
	try:
		saved = generate_heatmap.generate(args.database, heatmap_path)
		print(f"Heatmap saved: {saved}")
	except Exception as exc:
		print(f"Warning: could not generate heatmap: {exc}", file=sys.stderr)

	return 0


if __name__ == "__main__":
	raise SystemExit(main())
