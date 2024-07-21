import os
from Bio import SeqIO

AMINO_ACIDS = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}

class Fasta_Content:
  fasta_filename = str()
  """
  FASTA file name without directories and '.fasta'
  """
  fasta_content = dict()

  def __init__(self)->None:
    return

  def parse_fasta_file(self, file:str)->None:
    """
    Parses an opened FASTA file for all of its IDs and sequences
    """

    basename, extension = os.path.splitext(os.path.basename(file))
    
    # ensure that FASTA file has a FASTA extension
    if not extension.replace('.', '') in {"fasta", "fas", "fa", "fna", "ffn", "faa", "mpfa", "frn"}:
        raise ValueError

    # get FASTA file name
    self.fasta_filename = os.path.splitext(os.path.basename(file))[0]

    # parse FASTA file for sequence content
    for record in SeqIO.parse(file, 'fasta'):
      self.fasta_content[record.id] = record.seq

  def get_fasta_filename(self)->str:
    return self.fasta_filename

  def get_ids(self)->list:
    """
    Returns all IDs from a FASTA file.
    """
    return list(self.fasta_content.keys())

  def get_sequences(self)->list:
    """
    Returns all sequences from a FASTA file.
    """
    return list(self.fasta_content.values())

