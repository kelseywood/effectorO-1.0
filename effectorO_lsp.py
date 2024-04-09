from argparse import ArgumentParser
from subprocess import run

from Bio.Blast.Applications import NcbiblastpCommandline

def main():
  parser = ArgumentParser(prog="effectorO_lsp", description="Python script that run's the EffectorO LSP pipeline.")
  parser.add_argument("-d", "--ncbi_database", type=str, help="Input NCBI database", required=True)
  parser.add_argument("-i", "--input_fasta", type=str, help="Input FASTA file for queried genome", required=True)
  parser.add_argument("-n", "--genome_name", type=str, help="Genome name used for identification", required=True)
  args = parser.parse_args()

  database:str    = args.ncbi_database
  fasta:str       = args.input_fasta
  genome_name:str = args.genome_name

  NcbiblastpCommandline(cmd="blastp",
                        db=database,
                        query=fasta,
                        outfmt="6 std qcovs",
                        num_threads=4,
                        out=f"sp_{genome_name}_vs_all.tab")()


if __name__ == "__main__":
  main()