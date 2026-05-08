#!/usr/bin/env python
"""Round-trip verification for SDF files: encode to AMSR, decode, compute RMSD.

Each molecule is tested with a default encoding plus N_RANDOM_SEEDS randomized
encodings (same as the pytest suite).

Usage:
    python roundtrip_sdf.py ~/sdf              # sequential
    python roundtrip_sdf.py ~/sdf -j 10        # 10 parallel workers
    python roundtrip_sdf.py ~/sdf -j 10 -t 0.5 # stricter threshold
"""

import csv
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
from typing import Annotated, Optional

import typer

from amsr.roundtrip import SDF_RMSD_THRESHOLD, RoundtripSDF

app = typer.Typer(
    add_completion=False,
    context_settings={"help_option_names": ["-h", "--help"]},
)


def _sdf_work(input_dir: Path, n_seeds: int):
    """Yield (sdf_path, seed) lazily by recursively scanning directories."""
    seeds = [None] + list(range(n_seeds))
    for dirpath, dirnames, filenames in os.walk(input_dir):
        dirnames.sort()
        for fname in sorted(filenames):
            if fname.endswith(".sdf"):
                for seed in seeds:
                    yield os.path.join(dirpath, fname), seed


def _report(r, counts):
    counts[3] += 1
    if r["status"] == "PASSED":
        counts[0] += 1
        mark = "\033[32mPASSED\033[0m"
    elif r["status"] == "FAILED":
        counts[1] += 1
        mark = "\033[31mFAILED\033[0m"
    else:
        counts[2] += 1
        mark = "\033[33mERROR\033[0m"
    rmsd_str = f"rmsd={r['rmsd']:.3f}" if isinstance(r["rmsd"], float) else r.get("error", "")
    time_str = f"{r['time']:.2f}s" if r["time"] else ""
    n = counts[3]
    sys.stdout.write(f"[{n}] {r['name']}[seed{r['seed']}] {mark} {rmsd_str} {time_str}\n")
    if n % 100 == 0:
        sys.stdout.write(
            f"  --- {counts[0]} passed, {counts[1]} failed, {counts[2]} errors / {n} total ---\n"
        )
    sys.stdout.flush()


def _error_result(path, seed):
    return {
        "name": os.path.splitext(os.path.basename(path))[0],
        "seed": "0" if seed is None else str(seed + 1),
        "amsr": "",
        "rmsd": "",
        "status": "ERROR",
        "time": 0.0,
    }


@app.command()
def main(
    input_dir: Annotated[
        Path, typer.Argument(help="Directory (recursively searched) containing SDF files")
    ],
    output: Annotated[
        Optional[Path], typer.Option("-o", help="Output directory for SDF files (omit to skip)")
    ] = None,
    threshold: Annotated[float, typer.Option("-t", help="RMSD threshold")] = SDF_RMSD_THRESHOLD,
    jobs: Annotated[int, typer.Option("-j", help="Number of parallel workers")] = (
        os.cpu_count() or 1
    ),
    nseeds: Annotated[int, typer.Option("-n", help="Number of random seeds per molecule")] = 10,
):
    """Round-trip verification: encode each SDF to AMSR, decode, compute RMSD."""
    if not input_dir.is_dir():
        print(f"Error: {input_dir} is not a directory", file=sys.stderr)
        sys.exit(1)

    output_dir = str(output) if output is not None else None
    if output_dir is not None:
        os.makedirs(output_dir, exist_ok=True)

    csv_path = os.path.join(".", "roundtrip_sdf_out.csv")
    csv_file = open(csv_path, "w", newline="")
    csv_writer = csv.writer(csv_file)
    csv_writer.writerow(["name", "seed", "amsr", "rmsd", "status", "time_s"])

    counts = [0, 0, 0, 0]  # pass, fail, error, total

    def _record(r):
        _report(r, counts)
        rmsd_str = f"{r['rmsd']:.3f}" if isinstance(r["rmsd"], float) else ""
        csv_writer.writerow(
            [r["name"], r["seed"], r["amsr"], rmsd_str, r["status"], f"{r['time']:.3f}"]
        )
        csv_file.flush()

    if jobs <= 1:
        for path, seed in _sdf_work(input_dir, nseeds):
            _record(RoundtripSDF(path, seed, threshold, output_dir))
    else:
        max_pending = jobs * 2
        work = _sdf_work(input_dir, nseeds)
        with ProcessPoolExecutor(max_workers=jobs) as executor:
            futures = {}
            exhausted = False
            # Fill initial batch
            for path, seed in work:
                futures[executor.submit(RoundtripSDF, path, seed, threshold, output_dir)] = (
                    path,
                    seed,
                )
                if len(futures) >= max_pending:
                    break
            else:
                exhausted = True

            while futures:
                done = next(iter(as_completed(futures)))
                path, seed = futures.pop(done)
                try:
                    _record(done.result())
                except Exception:
                    _record(_error_result(path, seed))
                # Submit more work
                if not exhausted:
                    for path, seed in work:
                        futures[
                            executor.submit(RoundtripSDF, path, seed, threshold, output_dir)
                        ] = (path, seed)
                        if len(futures) >= max_pending:
                            break
                    else:
                        exhausted = True

    csv_file.close()
    total = counts[3]
    print(f"\n{'='*60}")
    print(f"{counts[0]} passed, {counts[1]} failed, {counts[2]} errors / {total} total")
    print(f"Results: {csv_path}")


if __name__ == "__main__":
    app()
