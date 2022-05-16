#!/usr/bin/env python

# Run e2 programs by running the commands with -h

from pathlib import Path
from multiprocessing import Pool
import subprocess
import sys


def run_prog(prog):
    proc = subprocess.run([prog, "-h"], stdout=subprocess.DEVNULL)
    print(f"Running: {' '.join(proc.args)}", flush=True)

    return prog if proc.returncode else None


def main():
    MYDIR = Path(__file__).parent
    PROGS_DIR = MYDIR.parent / "programs"

    progs = set(p.name for p in PROGS_DIR.glob('e2*.py'))

    with open(MYDIR / "programs_no_test.txt", "r") as fin:
        progs_exclude = set()
        for line in fin:
            progs_exclude.add(line.split()[0])

    print("\nRemoving programs from test list...")
    for f in progs_exclude:
        print(f"... {f}")

    progs -= progs_exclude

    with Pool() as pool:
        failed_progs = [p for p in pool.map(run_prog, progs) if p]

    print(f"\nTotal failed programs: {len(failed_progs)} / {len(progs)}")
    for prog in failed_progs:
        print(prog)
    print()

    if failed_progs:
        sys.exit(1)


if __name__ == '__main__':
    main()
