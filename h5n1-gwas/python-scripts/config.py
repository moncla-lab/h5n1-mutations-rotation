tree_path = "/Users/jort/coding/h5n1-mutations-rotation/h5n1-gwas/test-data/flu_avian_h5n1_pb2.json"
#tree_path = "/Users/jort/coding/h5n1-mutations-rotation/base-build/auspice/flu_avian_h5n1_pb2.json"

"""For the tree you are reading in, you need to specify which attribute encodes the host value. For the avian 
flu trees on nextstrain, this attribute is `host`. If you specified a different label, like `host_species` or 
`Host`, you would need to change this. You can check this by manually looking at the tree json file in a text
editor, or by parsing through a tree and reading the attributes. """
host_annotation = 'host'

"""specify the gene you are running, which hosts to compare, method, and minimum required count. Here, host 1 is 
the host that we want to find mutations that are enriched in, host 2 is the background host"""
gene = "PB2"
host1 = "Human"
host2 = "Avian"
minimum_required_count = 0
method = "counts"


## total iterations = iterations * # CPU cores (retrieved by mp.cpu_count() – 10 on Jordan's laptop)
iterations = 1000