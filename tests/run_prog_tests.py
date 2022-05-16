#!/usr/bin/env python

# Run e2 programs by running the commands with -h

from pathlib import Path
import subprocess
import sys


failed_progs = []


def run_prog(prog):
    proc = subprocess.run([prog, "-h"], stdout=subprocess.DEVNULL)
    print(f"Running: {' '.join(proc.args)}", flush=True)
    if proc.returncode:
        failed_progs.append(prog)


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

    for prog in progs:
        run_prog(prog)

    print(f"\nTotal failed programs: {len(failed_progs)} / {len(progs)}")
    for prog in failed_progs:
        print(prog)

    if failed_progs:
        sys.exit(1)


if __name__ == '__main__':
    main()
