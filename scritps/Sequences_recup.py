from Bio import Entrez, SeqIO
from tqdm import tqdm
import os
import time

def fetch_ncbi_genomes_by_genus(
    genus,
    outdir,
    email,
    retmax=20,
    delay=0.4
):
    """
    Download complete bacterial genomes from NCBI by genus.
    """

    Entrez.email = email
    Entrez.tool = "comparative_genomics_demo"

    os.makedirs(outdir, exist_ok=True)

    query = f"{genus}[Organism] AND complete genome[Title]"

    print(f"Searching NCBI for: {query}")

    search_handle = Entrez.esearch(
        db="nucleotide",
        term=query,
        retmax=retmax
    )
    search_results = Entrez.read(search_handle)
    search_handle.close()

    ids = search_results["IdList"]

    print(f"{len(ids)} genomes found")

    for seq_id in tqdm(ids, desc="Downloading genomes"):
        try:
            fetch_handle = Entrez.efetch(
                db="nucleotide",
                id=seq_id,
                rettype="fasta",
                retmode="text"
            )

            record = SeqIO.read(fetch_handle, "fasta")
            fetch_handle.close()

            outfile = os.path.join(outdir, f"{record.id}.fasta")
            SeqIO.write(record, outfile, "fasta")

            time.sleep(delay)  # respect NCBI servers

        except Exception as e:
            print(f"Failed to fetch {seq_id}: {e}")

if __name__ == "__main__":

    outdir = "data/genomes/escherichia"

    fetch_ncbi_genomes_by_genus(
        genus="Escherichia",
        outdir=outdir,
        email="your.email@domain.com",
        retmax=20
    )
