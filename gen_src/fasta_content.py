from argparse import FileType
from Bio.SeqIO import parse

AMINO_ACIDS = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}

class Fasta_Content:
  fasta_filename = str()
  """
  FASTA file name without directories and '.fasta'
  """
  fasta_content = dict()

  def __init__(self)->None:
    return

  def parse_fasta_file(self, file:FileType)->None:
    """
    Parses an opened FASTA file for all of its IDs and sequences
    """

    # get FASTA file name
    fasta_filename = str(file.name)
    self.fasta_filename = fasta_filename[fasta_filename.rfind('/')+1:fasta_filename.find('.')]

    # parse FASTA file for sequence content
    for record in parse(file, 'fasta'):
      self.fasta_content[record.id] = record.seq

  def get_fasta_filename(self)->str:
    return self.fasta_filename

  def get_ids(self)->list:
    """
    Returns all IDs from a FASTA file.
    """
    return self.fasta_content.keys()

  def get_sequences(self)->list:
    """
    Returns all sequences from a FASTA file.
    """
    return self.fasta_content.values()
