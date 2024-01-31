from argparse import FileType

AMINO_ACIDS = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'}

class Fasta_Content:
  fasta_filename = str()
  """
  FASTA file name without directories and '.fasta'
  """
  fasta_content = list()
  """
  A list of 3-item tuples that contain
  (ID, header content list (not including ID), sequence)
  """

  def __init__(self)->None:
    return
  
  def parse_fasta_file(self, file:FileType('r'))->None:
    """
    Parses an opened FASTA file for all of its IDs and sequences
    """
    
    # get FASTA file name
    fasta_filename = str(file.name)
    self.fasta_filename = fasta_filename[fasta_filename.rfind('/')+1:fasta_filename.find('.')]
    
    # parse FASTA file for sequence content
    while True:
      line = next(file, None)
      if line == None: break

      if line[0] == '>':
        # get ID and header
        id = line[1:].split()[0]
        header = line[len(id)+1:].strip().split()

        # perform actions on sequence before import
        sequence = next(file, None)
        if sequence == None: raise ValueError("An error occurred: A sequence was not associated with an ID")

        with open("warnings.log", 'w') as fwarn:
          for i in reversed(range(len(sequence))):
            if not sequence[i].upper() in AMINO_ACIDS:
              fwarn.write(f"Warning: character {sequence[i]} found in location {i} of sequence {id}. Deleting...\n")
              sequence = sequence[:i] + sequence[i+1:]

        # add sequence info to class
        self.fasta_content.append((id, header, sequence))
  
  def get_fasta_filename(self)->str:
    return self.fasta_filename

  def get_ids(self)->list:
    """
    Returns all IDs from a FASTA file.
    """
    return [id for id, _, _ in self.fasta_content]

  def get_headers(self)->list:
    """
    Returns all header contents from a FASTA file.
    """
    return [header for _, header, _ in self.fasta_content]

  def get_sequences(self)->list:
    """
    Returns all sequences from a FASTA file.
    """
    return [sequence for _, _, sequence in self.fasta_content]
