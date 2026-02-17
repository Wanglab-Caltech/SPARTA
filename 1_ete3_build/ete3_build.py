#!/usr/bin/env python3
"""
Create an ETE3 NCBI taxonomy SQLite database at taxdump/ncbi.sqlite.

Usage:
  python ete3_build.py
"""

from __future__ import annotations
import argparse
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser(description="Download/build ETE3 NCBI taxonomy SQLite database.")
    parser.add_argument("--outdir", default="taxdump", help="Output directory (default: taxdump)")
    parser.add_argument("--dbfile", default="ncbi.sqlite", help="SQLite filename (default: ncbi.sqlite)")
    args = parser.parse_args()

    outdir = Path(args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    dbpath = outdir / args.dbfile

    try:
        from ete3 import NCBITaxa
    except ImportError as e:
        raise SystemExit(
            "ete3 is not installed in this environment.\n"
            "Install with: conda install -c conda-forge ete3\n"
            f"Original error: {e}"
        )

    print(f"Building/updating ETE3 taxonomy DB at: {dbpath}")
    ncbi = NCBITaxa(dbfile=str(dbpath))
    ncbi.update_taxonomy_database()
    print("Done:", ncbi.dbfile)


if __name__ == "__main__":
    main()
