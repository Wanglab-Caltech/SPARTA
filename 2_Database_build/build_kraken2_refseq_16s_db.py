#!/usr/bin/env python3
"""
Build a custom Kraken2 16S database from NCBI RefSeq TargetedLoci 16S FASTAs
(Bacteria + Archaea). TaxIDs are attached using nucl_gb.accession2taxid.gz,
producing headers like: >NR_XXXX|kraken:taxid|12345
"""

from __future__ import annotations
import argparse
import gzip
import shutil
import subprocess
from pathlib import Path


ARCHAEA_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Archaea/archaea.16SrRNA.fna.gz"
BACTERIA_URL = "https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz"
ACC2TAXID_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz"


def run(cmd: list[str], cwd: Path | None = None) -> None:
    print("\n>>> " + " ".join(cmd))
    subprocess.run(cmd, cwd=str(cwd) if cwd else None, check=True)


def ensure_exe(name: str) -> None:
    if shutil.which(name) is None:
        raise SystemExit(f"Missing executable in PATH: {name}")


def accession_from_fasta_header(line: str) -> str:
    tok = line[1:].strip().split()[0]
    if "." in tok:
        tok = tok.split(".", 1)[0]
    return tok


def write_combined_fasta(out_fasta: Path, gz_files: list[Path]) -> None:
    with out_fasta.open("wt") as out:
        for gz in gz_files:
            with gzip.open(gz, "rt") as fh:
                shutil.copyfileobj(fh, out)


def extract_accessions(fasta_path: Path, out_txt: Path) -> list[str]:
    accs = set()
    with fasta_path.open("rt") as fh:
        for line in fh:
            if line.startswith(">"):
                accs.add(accession_from_fasta_header(line))
    acc_list = sorted(accs)
    out_txt.write_text("\n".join(acc_list) + "\n")
    return acc_list


def build_ref2taxid(accessions: set[str], acc2taxid_gz: Path, out_tsv: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with gzip.open(acc2taxid_gz, "rt") as fh, out_tsv.open("wt") as out:
        _ = next(fh, None)  # skip header line
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            acc = parts[0]
            if acc in accessions:
                taxid = parts[2]
                mapping[acc] = taxid
                out.write(f"{acc}\t{taxid}\n")
    return mapping


def write_fixed_headers(in_fasta: Path, out_fasta: Path, acc2taxid: dict[str, str]) -> None:
    with in_fasta.open("rt") as inp, out_fasta.open("wt") as out:
        for line in inp:
            if line.startswith(">"):
                acc = accession_from_fasta_header(line)
                taxid = acc2taxid.get(acc)
                if taxid:
                    out.write(f">{acc}|kraken:taxid|{taxid}\n")
                else:
                    out.write(line)
            else:
                out.write(line)


def cleanup_taxonomy(dbdir: Path, keep_extra: bool) -> None:
    taxdir = dbdir / "taxonomy"
    if not taxdir.exists():
        print(f"WARNING: taxonomy dir not found: {taxdir}")
        return

    keep = ["prelim_map.txt", "names.dmp", "nodes.dmp"]
    if keep_extra:
        keep += ["merged.dmp", "delnodes.dmp"]

    existing = {p.name for p in taxdir.iterdir() if p.is_file()}
    keep_present = [k for k in keep if k in existing]

    print("\nCleaning taxonomy folder. Keeping:")
    for k in keep_present:
        print(" -", k)

    tmp = taxdir / "_keep_tmp"
    tmp.mkdir(exist_ok=True)

    for k in keep_present:
        shutil.move(str(taxdir / k), str(tmp / k))

    for item in taxdir.iterdir():
        if item.name == "_keep_tmp":
            continue
        if item.is_dir():
            shutil.rmtree(item)
        else:
            item.unlink(missing_ok=True)

    for item in tmp.iterdir():
        shutil.move(str(item), str(taxdir / item.name))
    tmp.rmdir()

    print("\nFinal taxonomy directory contents:")
    run(["ls", "-lh", str(taxdir)])


def build_bracken_indices(dbdir: Path, threads: int, kmer: int, readlen: int) -> None:
    """
    Run bracken-build and place outputs into:
      <dbdir>/bracken_indices_<readlen>/
        database<readlen>mers.kmer_distrib
        database<readlen>mers.kraken
    """
    run(["bracken-build", "-d", str(dbdir), "-t", str(threads), "-k", str(kmer), "-l", str(readlen)])

    src_kmer = dbdir / f"database{readlen}mers.kmer_distrib"
    src_kra = dbdir / f"database{readlen}mers.kraken"

    outdir = dbdir / f"bracken_indices_{readlen}"
    outdir.mkdir(parents=True, exist_ok=True)

    if not src_kmer.exists():
        raise SystemExit(f"Expected Bracken output not found: {src_kmer}")
    shutil.move(str(src_kmer), str(outdir / src_kmer.name))

    if src_kra.exists():
        shutil.move(str(src_kra), str(outdir / src_kra.name))
    else:
        print(f"NOTE: {src_kra.name} not found; continuing (some Bracken builds only create .kmer_distrib)")

    print("\nBracken indices directory contents:")
    run(["ls", "-lh", str(outdir)])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--db", default="16s_db", help="Kraken2 DB directory (default: 16s_db)")
    parser.add_argument("--threads", type=int, default=8, help="Threads (default: 8)")
    parser.add_argument("--workdir", default="build_16s_refseq", help="Work directory (default: build_16s_refseq)")
    parser.add_argument("--no-masking", action="store_true", default=True, help="Use --no-masking (default: true)")
    parser.add_argument(
        "--keep-taxonomy-extra",
        action="store_true",
        help="Also keep merged.dmp and delnodes.dmp (recommended). Default keeps only prelim_map/names/nodes.",
    )
    parser.add_argument("--bracken-k", type=int, default=35, help="Bracken k-mer length (default: 35)")
    # indices can be adjusted based on the expected read length of the sequencing
    parser.add_argument("--bracken-l", type=int, default=1200, help="Bracken read length (default: 1200)")
    parser.add_argument("--skip-bracken", action="store_true", help="Skip Bracken index build step")
    args = parser.parse_args()

    ensure_exe("kraken2-build")
    ensure_exe("wget")
    if not args.skip_bracken:
        ensure_exe("bracken-build")

    workdir = Path(args.workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)

    dbdir = Path(args.db).resolve()
    print(f"\nWorking directory: {workdir}")
    print(f"Kraken2 DB directory: {dbdir}")

    # 1) download RefSeq TargetedLoci 16S FASTAs
    run(["wget", "-c", ARCHAEA_URL], cwd=workdir)
    run(["wget", "-c", BACTERIA_URL], cwd=workdir)
    archaea_gz = workdir / "archaea.16SrRNA.fna.gz"
    bacteria_gz = workdir / "bacteria.16SrRNA.fna.gz"

    combined_fasta = workdir / "16S.fasta"
    print(f"\nCombining FASTAs -> {combined_fasta}")
    write_combined_fasta(combined_fasta, [archaea_gz, bacteria_gz])

    # 2) extract accession IDs
    accessions_txt = workdir / "accessions.txt"
    acc_list = extract_accessions(combined_fasta, accessions_txt)
    accessions = set(acc_list)
    print(f"Extracted {len(accessions)} unique accessions -> {accessions_txt}")

    # 3) download accession->taxid mapping
    run(["wget", "-c", ACC2TAXID_URL], cwd=workdir)
    acc2taxid_gz = workdir / "nucl_gb.accession2taxid.gz"

    ref2taxid_tsv = workdir / "ref2taxid.tsv"
    print(f"\nFiltering accession2taxid -> {ref2taxid_tsv} (streaming; may take a bit)")
    acc2taxid = build_ref2taxid(accessions, acc2taxid_gz, ref2taxid_tsv)
    print(f"Mapped {len(acc2taxid)} accessions to taxids")

    # 4) rewrite FASTA headers
    fixed_headers = workdir / "fixed_headers.fasta"
    print(f"\nWriting fixed headers -> {fixed_headers}")
    write_fixed_headers(combined_fasta, fixed_headers, acc2taxid)

    # 5) build Kraken2 DB
    if dbdir.exists():
        print(f"\nRemoving existing DB: {dbdir}")
        shutil.rmtree(dbdir)
    dbdir.mkdir(parents=True, exist_ok=True)

    run(
        [
            "kraken2-build",
            "--download-taxonomy",
            "--use-ftp",
            "--skip-maps",
            "--db",
            str(dbdir),
            "--threads",
            str(args.threads),
        ]
    )

    add_cmd = ["kraken2-build", "--add-to-library", str(fixed_headers), "--db", str(dbdir)]
    if args.no_masking:
        add_cmd.append("--no-masking")
    run(add_cmd)

    build_cmd = ["kraken2-build", "--build", "--db", str(dbdir), "--threads", str(args.threads)]
    if args.no_masking:
        build_cmd.append("--no-masking")
    run(build_cmd)

    # 6) sanity check
    run(["kraken2-inspect", "--db", str(dbdir), "--skip-counts"])

    # 7) Ensure names/nodes exist
    nodes = dbdir / "taxonomy" / "nodes.dmp"
    names = dbdir / "taxonomy" / "names.dmp"
    if not nodes.exists() or not names.exists():
        raise SystemExit(f"Taxonomy files missing: {nodes} or {names}")
    run(["ls", "-lh", str(nodes), str(names)])

    # 8) garbage cleanup
    cleanup_taxonomy(dbdir, keep_extra=args.keep_taxonomy_extra)

    # 9) Bracken build + move into bracken_indices_<L>/
    if not args.skip_bracken:
        build_bracken_indices(dbdir, threads=args.threads, kmer=args.bracken_k, readlen=args.bracken_l)

    print("\n✅ Done.")
    print(f"DB built at: {dbdir}")
    print(f"Workdir kept at: {workdir}")


if __name__ == "__main__":
    main()
