from Bio import SeqIO
import pandas as pd
import os

def compute_gc(seq):
    seq = seq.upper()
    gc = seq.count("G") + seq.count("C")
    return (gc / len(seq)) * 100 if len(seq) > 0 else 0


def extract_genome_features(genome_dir, genus, output_csv):
    records = []

    for file in os.listdir(genome_dir):
        if not file.endswith(".fasta"):
            continue

        filepath = os.path.join(genome_dir, file)
        contigs = list(SeqIO.parse(filepath, "fasta"))

        genome_size = sum(len(rec.seq) for rec in contigs)
        full_sequence = "".join(str(rec.seq) for rec in contigs)
        gc_content = compute_gc(full_sequence)

        records.append({
            "genome_id": file.replace(".fasta", ""),
            "genus": genus,
            "genome_size_bp": genome_size,
            "gc_content_percent": round(gc_content, 2),
            "n_contigs": len(contigs)
        })

    df = pd.DataFrame(records)
    df.to_csv(output_csv, index=False)

    print(f"Saved {len(df)} genomes to {output_csv}")

if __name__ == "__main__":

    genome_dir = "data/genomes/escherichia"
    output_csv = "results/genome_features.csv"

    extract_genome_features(
        genome_dir=genome_dir,
        genus="Escherichia",
        output_csv=output_csv
    )
