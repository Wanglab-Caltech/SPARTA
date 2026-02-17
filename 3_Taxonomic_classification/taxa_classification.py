"""
This script takes a folder of FASTA/FASTQ files, runs Kraken2 on each one to do taxonomic classification, and saves the Kraken report.
It then runs Bracken on those Kraken reports to estimate abundances (at the level you set, like species),
and converts the Bracken results into a clean CSV that includes the full taxonomy lineup (Kingdom → Species) using an ETE3 NCBI taxonomy database.
"""
import os
import sys
import subprocess
import glob
import pandas as pd
from ete3 import NCBITaxa

def run_kraken2_and_bracken(input_dir, kraken_db='16s_db', bracken_db='16s_db/bracken_indices_1200', taxdump_db='taxdump/ncbi.sqlite', min_abundance=0.0001):
    input_dir = os.path.abspath(input_dir)
    output_dir = os.path.join(input_dir, "kraken2_results")
    os.makedirs(output_dir, exist_ok=True)

    # Accept FASTA and FASTQ (gzipped or not)
    exts = (".fastq", ".fastq.gz", ".fq", ".fq.gz", ".fasta", ".fasta.gz", ".fa", ".fa.gz", ".fna", ".fna.gz")
    seq_files = [f for f in os.listdir(input_dir) if f.endswith(exts)]

    if not seq_files:
        print("No FASTA/FASTQ files found in:", input_dir)
        return

    print(f"Found {len(seq_files)} sequence file(s). Running Kraken2...\n")

    for seq_file in seq_files:
        input_path = os.path.join(input_dir, seq_file)
        base_name = os.path.splitext(os.path.splitext(seq_file)[0])[0]

        output_tsv = os.path.join(output_dir, f"kraken_output_{base_name}.tsv")
        output_report = os.path.join(output_dir, f"kraken_report_{base_name}.txt")

        cmd = [
            "kraken2",
            "--db", kraken_db,
            "--use-names",
            "--report", output_report,
            "--output", output_tsv,
            input_path
        ]

        print(f">>> Running Kraken2 on {seq_file}")
        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            print(f"Error processing {seq_file}: {e}")
        print(f"✓ Output saved to: {output_tsv}, {output_report}\n")

        # --- NEW: delete kraken output tsv (not used downstream) ---
        try:
            os.remove(output_tsv)
        except FileNotFoundError:
            pass

    print("✅ Kraken2 analysis completed.\n")

    print("Running Bracken on each kraken_report_*.txt...\n")

    for file in os.listdir(output_dir):
        # --- CHANGED: skip bracken-generated report copies ---
        if file.startswith("kraken_report_") and file.endswith(".txt") and "_bracken_" not in file:
            report_path = os.path.join(output_dir, file)
            base_name = os.path.splitext(file)[0].replace("kraken_report_", "")
            bracken_output = os.path.join(output_dir, f"bracken_output_{base_name}.tsv")

            bracken_cmd = [
                "bracken",
                "-d", bracken_db,
                "-i", report_path,
                "-o", bracken_output,
                "-r", "1200",
                "-l", "S",
                "-t", "2"
            ]

            print(f">>> Running Bracken on {file}")
            try:
                subprocess.run(bracken_cmd, check=True)
            except subprocess.CalledProcessError as e:
                print(f"Error running Bracken on {file}: {e}")
            print(f"✓ Bracken output: {bracken_output}\n")

            # --- NEW: delete bracken-generated report copies like *_bracken_species.txt ---
            for junk in glob.glob(os.path.join(output_dir, f"kraken_report_{base_name}_bracken_*.txt")):
                try:
                    os.remove(junk)
                except FileNotFoundError:
                    pass

    print("✅ Bracken analysis completed for all reports.\n")

    print("Converting Bracken TSV to CSV with taxonomic info...\n")

    ncbi = NCBITaxa(dbfile=taxdump_db)
    ranks_wanted = ["superkingdom", "phylum", "class", "order", "family", "genus", "species"]
    rank_to_column = {
        "superkingdom": "Kingdom",
        "phylum": "Phylum",
        "class": "Class",
        "order": "Order",
        "family": "Family",
        "genus": "Genus",
        "species": "Species"
    }

    converted_dir = os.path.join(output_dir, "bracken_converted_csv")
    os.makedirs(converted_dir, exist_ok=True)
    tsv_files = sorted(
        p for p in glob.glob(os.path.join(output_dir, "bracken_output_*.tsv"))
        if "_bracken_" not in os.path.basename(p)
    )

    print(f"Found {len(tsv_files)} Bracken TSV files")

    for tsv_path in tsv_files:
        base = os.path.basename(tsv_path).replace(".tsv", "")
        out_csv_path = os.path.join(converted_dir, f"{base}.csv")

        df = pd.read_csv(tsv_path, sep="\t")
        df = df[df["fraction_total_reads"] >= min_abundance].copy()

        records = []
        for _, row in df.iterrows():
            taxid = int(row["taxonomy_id"])
            abundance = row["fraction_total_reads"]
            try:
                lineage = ncbi.get_lineage(taxid)
                names = ncbi.get_taxid_translator(lineage)
                ranks = ncbi.get_rank(lineage)

                record = {
                    "Taxonomy ID": taxid,
                    "Relative Abundance": abundance
                }
                for tid in lineage:
                    rank = ranks.get(tid)
                    if rank in ranks_wanted:
                        record[rank_to_column[rank]] = names[tid]
                records.append(record)
            except Exception:
                print(f"⚠️ TaxID {taxid} could not be resolved")

        df_out = pd.DataFrame(records)
        for col in ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]:
            if col not in df_out.columns:
                df_out[col] = ""

        df_out = df_out[["Taxonomy ID", "Relative Abundance", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]]
        df_out.to_csv(out_csv_path, index=False)
        print(f"✅ Converted and saved: {out_csv_path}")

    print("\nAll processing complete.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python batch_kraken2_bracken_converter.py <seq_directory>")
        sys.exit(1)

    run_kraken2_and_bracken(sys.argv[1])
