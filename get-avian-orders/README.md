# README

FASTA ID formatting:

{strain}|{virus}|{isolate_id}|{date}|{region}|{country}|{division}|{location}|{host}|{subtype}|{originating_lab}|{submitting_lab}|{h5_clade}<br/><br/>

This python script parses FASTA files, scanning for IDs with {host} == 'avian', identifying their avian {host_species} as indicated in the {strain} name [formatted as {type}/{host_species}/{location}/{isolate}/{year}], and determining a standardized species name based on avian-species-synonyms.txt. The taxonomical order of that standardized species is then determined with avian-species-orders.txt (currently only contains labels 'a' for anseriformes, 'g' for galliformes, 'o' for other). The script finally will write a new FASTA file with IDs in which {host} is replaced with the order name for any sequences that could be classified. Any IDs that could not be classified (non-avian host, missing in txt files, or labeled as 'o') will not be altered; this file will be named as {original_file_name}_avian_orders.fasta.

To use the script, change <b>"/Users/jort/coding/AvianOrders/"</b> to the directory containing your FASTA file(s) and <b>["h5n1_ha.fasta", "h5n1_pb2.fasta"]</b> to a list containing the names of your file(s).
