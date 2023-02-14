from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd


## specify directory containing fasta files
## output files will be saved here (with '_clades' appended to file name)
fasta_dir = '/Users/jort/Desktop/test/'

## specify names of fasta files
fasta_file_names = ['h5n1_ha.fasta', 'h5n1_pb2.fasta']

## specify path to h5nx-clades.tsv
tsv_path = '/Users/jort/coding/h5n1-mutations-rotation/get-h5n1-clade/h5nx-clades.tsv'


##### user input above #####


## make dataframe from clade tsv file
clade_df = pd.read_table(tsv_path)

## functions for re-writing FASTA file
def get_clade_id(name):
  '''get FASTA ID with H5N1 clade'''
  clade = ""
  name_comps = name.split("|")
  if name_comps[0] in clade_df.name.values:
    clade = str(clade_df[clade_df['name'] == name_comps[0]].clade.tolist()[0])
  if clade != "":
    name_comps[12] = clade
  return("|".join(name_comps))

def write_order_fasta(fasta_file_name):
  '''write FASTA file with new IDs'''
  new_file_name = "".join([fasta_file_name.split(".")[0], "_clades.", fasta_file_name.split(".")[1]])
  fasta_sequences = SeqIO.parse(open(fasta_file_name),'fasta')
  records = []
  for fasta in fasta_sequences:
    name, sequence = fasta.id, fasta.seq
    new_name = get_clade_id(name)
    record = SeqRecord(sequence, new_name, "", "")
    records.append(record)
  SeqIO.write(records, new_file_name, "fasta")

## write new FASTA file for all files in list
for file in fasta_file_names:
  file_path = fasta_dir+file
  write_order_fasta(file_path)