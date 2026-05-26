#!/usr/bin/env python
# e3ctomo_gain.py  --  Continuous-tilt gain reference utilities.
from past.utils import old_div
from builtins import range
from EMAN2 import *
import sys, os, shutil

def main():
    progname = os.path.basename(sys.argv[0])
    usage = """e3ctomo_gain.py [options]

Continuous-tilt gain reference utilities.

  importgain : copy/convert an existing gain reference into the project.
  creategain : compute a new gain reference from a set of movie files.
"""
    parser = EMArgumentParser(usage=usage, version=EMANVERSION)

    # importgain mode
    parser.add_argument("--gainsrc", type=str, default="", help="Input gain reference file to import.", guitype='filebox', browser="EMBrowserWidget(withmodal=True, multiselect=False)", filecheck=False, row=0, col=0, rowspan=1, colspan=3, mode="importgain")

    # creategain mode
    parser.add_pos_argument(name="movies", help="Movie files to compute the gain reference from.", default="", guitype='filebox', browser="EMBrowserWidget(withmodal=True, multiselect=True)", row=0, col=0, rowspan=1, colspan=3, mode="creategain")
    parser.add_argument("--gainname", type=str, default="GainFiducial.hdf", help="Output filename for the computed gain reference.", guitype='strbox', row=2, col=0, rowspan=1, colspan=3, mode="creategain")

    # shared: importgain and creategain
    parser.add_argument("--dest", type=str, default=".", help="Destination directory.", guitype='filebox', browser="EMBrowserWidget(withmodal=True, multiselect=False)", filecheck=False, row=1, col=0, rowspan=1, colspan=3, mode="importgain,creategain")
    parser.add_header(name="apixhdr", help="", title=u"Å/pixel", row=3, col=0, rowspan=1, colspan=1, mode="importgain,creategain")
    parser.add_argument("--apix", type=float, default=1.0, help=u"Å/pixel — Angstroms per pixel.", guitype='floatbox', row=3, col=1, rowspan=1, colspan=1, mode="importgain,creategain")

    parser.add_argument("--ppid", type=int, default=-1, help="Set the PID of the parent process, used for cross-platform PPID.")

    (options, args) = parser.parse_args()

    logid = E2init(sys.argv, options.ppid)

    if options.gainsrc:
        # importgain: run e2proc2d.py to copy and set apix
        dest_dir = options.dest if options.dest else "."
        os.makedirs(dest_dir, exist_ok=True)
        dst = os.path.join(dest_dir, os.path.basename(options.gainsrc))
        cmd = f"e2proc2d.py {options.gainsrc} {dst} --apix={options.apix} --compressbits=4"
        print(cmd)
        ret = os.system(cmd)
        if ret != 0:
            print(f"e2proc2d.py exited with code {ret}")
    elif args:
        # creategain: run e3movie.py on supplied movies
        dest_dir = options.dest if options.dest else "."
        os.makedirs(dest_dir, exist_ok=True)
        gainname = options.gainname if options.gainname else "GainFiducial.hdf"
        est_gain = os.path.join(dest_dir, gainname)
        cmd = "e3movie.py " + " ".join(args) + f" --est_gain={est_gain}"
        print(cmd)
        ret = os.system(cmd)
        if ret != 0:
            print(f"e3movie.py exited with code {ret}")

    E2end(logid)

if __name__ == "__main__":
    main()
