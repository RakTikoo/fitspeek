"""Command-line interface for fitspeek."""

from __future__ import annotations

import argparse
import sys

from .core import make


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="fitspeek",
        description="Turn a pile of FITS image files into contact-sheet PNGs "
        "with header stamps, for quick triage.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    make_parser = sub.add_parser(
        "make", help="build contact sheets from FITS files or folders"
    )
    make_parser.add_argument(
        "inputs",
        nargs="+",
        help="FITS files and/or folders containing FITS files",
    )
    make_parser.add_argument(
        "--out",
        default="sheets",
        help="output folder for sheet PNGs and index.csv (default: sheets)",
    )
    make_parser.add_argument(
        "--cols",
        type=int,
        default=3,
        help="thumbnails per row and column on each sheet (default: 3)",
    )

    args = parser.parse_args(argv)
    if args.command == "make":
        if args.cols < 1 or args.cols > 10:
            parser.error("--cols must be in the range 1 to 10")
        make(args.inputs, args.out, args.cols)
        return 0
    return 2


if __name__ == "__main__":
    sys.exit(main())
