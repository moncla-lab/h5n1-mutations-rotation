from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os

## change to directory with fasta files
os.chdir("/Users/jort/coding/AvianOrders/")

## specify names of fasta files
fasta_file_names = ["h5n1_ha.fasta", "h5n1_pb2.fasta"]



## make dictionary to standardize species name
species_dict = {}

with open('avian-species-synonyms.txt', 'r') as synfile:
    lines = synfile.read()

for entry in lines.split("\n"):
  entry_list = entry.split("\t")
  if(len(entry_list)) > 1 and entry_list[0][0] != "#":
    species_dict[entry_list[0]] = entry_list[1]


## make dictionary to get avian order
orders_dict = {}

with open('avian-species-orders.txt', 'r') as orderfile:
    lines = orderfile.read()

for entry in lines.split("\n"):
  entry_list = entry.split("\t")
  if entry_list[1] == "a":
    orders_dict[entry_list[0]] = "anseriforme"
  elif entry_list[1] == "g":
    orders_dict[entry_list[0]] = "galliforme"
  elif entry_list[1] == "o":
    orders_dict[entry_list[0]] = "avian"


## get FASTA id with avian order as host
def get_order_name(name):
  order = ""
  name_comps = name.split("|")
  if name_comps[8] == "avian":
    nonstd_name = name.split("|")[0].split("/")[1].lower()
    if nonstd_name in species_dict:
      std_name = species_dict[nonstd_name]
      if std_name in orders_dict:
        order = orders_dict[std_name]
  if order != "":
    name_comps[8] = order
  return("|".join(name_comps))


## write FASTA file with new ids
def write_order_fasta(fasta_file_name):
  new_file_name = "".join([fasta_file_name.split(".")[0], "_avian_orders.", fasta_file_name.split(".")[1]])
  fasta_sequences = SeqIO.parse(open(fasta_file_name),'fasta')
  records = []
  for fasta in fasta_sequences:
    name, sequence = fasta.id, fasta.seq
    new_name = get_order_name(name)
    record = SeqRecord(sequence, new_name, "", "")
    records.append(record)
  SeqIO.write(records, new_file_name, "fasta")


## write new FASTA file for all files in list
for name in fasta_file_names:
  write_order_fasta(name)