## Specify the local source file containing your Auspice JSON tree
#tree_path = "/Users/jort/coding/h5n1-mutations-rotation/h5n1-gwas/test-data/flu_avian_h5n1_pb2.json"
tree_path = "/Users/jort/coding/h5n1-mutations-rotation/base-build/auspice/flu_avian_h5n1_pb2.json"
#tree_path = "/Users/jort/coding/h5n1-mutations-rotation/base-build/auspice/flu_avian_h5n1_ha.json"


## There are a couple different tools that can be used to parse trees, but my preferred is baltic. Baltic is a tool written
## in python by Gytis Dudas, and available here: https://github.com/evogytis/baltic. If installing via pip (recommended), import
## in tree_manager.py with `import baltic as bt`. Otherwise, import with imp from a local source file specified in the config file.
##
## If baltic is installed with pip, set baltic_path to "pip"
## If baltic is installed locally, set baltic_path to the local source file
baltic_path = "/Users/jort/coding/baltic/baltic/baltic.py" # or "pip"

## For the tree you are reading in, you need to specify which attribute encodes the host value. For the avian 
## flu trees on nextstrain, this attribute is `host`. If you specified a different label, like `host_species` or 
## `Host`, you would need to change this. You can check this by manually looking at the tree json file in a text
## editor, or by parsing through a tree and reading the attributes.
host_annotation = "host"

## Specify the gene you are running, which hosts to compare, method, and minimum required count. Here, host 1 is 
## the host that we want to find mutations that are enriched in, host 2 is the background host
gene = "PB2"
host1 = "Human"
host2 = "Avian"
method = "counts"
minimum_required_count = 0

## Specify the number of simulations to perform
iterations = 200000