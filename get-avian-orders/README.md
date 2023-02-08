# README

FASTA ID formatting:

<b>{strain}</b>|{virus}|{isolate_id}|{date}|{region}|{country}|{division}|{location}|<b>{host}</b>|{subtype}|{originating_lab}|{submitting_lab}|{h5_clade}<br/><br/>

This python script parses FASTA files, scanning for IDs with <b>{host}</b> == 'avian', identifying their avian <b>{host_species}</b> as indicated in the <b>{strain}</b> name [formatted as {type}/<b>{host_species}</b>/{location}/{isolate}/{year}], and determining a standardized species name based on avian-species-synonyms.txt. The taxonomical order of that standardized species is then determined with avian-species-orders.txt (currently only contains labels 'a' for anseriformes, 'g' for galliformes, 'o' for other). The script finally will write a new FASTA file with IDs in which <b>{host}</b> is replaced with the order name for any sequences that could be classified; this file will be named as OriginalFileName_avian_orders.fasta. Any IDs that could not be classified (non-avian host, missing in txt files, or labeled as 'o') will not be altered.

To use the script, change <b>"/Users/jort/coding/AvianOrders/"</b> to the directory containing your FASTA file(s) and <b>["h5n1_ha.fasta", "h5n1_pb2.fasta"]</b> to a list containing the names of your file(s).
