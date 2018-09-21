I am making the rplB tree for use in the circular figure. These proteins are taken directly from NCBI and
are present in two fasta files. The two files are full annotation and renamed. because the muscle output
can't deal with long fasta headers, I have shortened each fasta header to just the genus_species. The 
full annotation file contains the actual protein accession, etc.

using the fasta file, I created an alignment with the web tool MUSCLE in this dir. I used the file rplB_renamed.faa, which is a fasta
file of the rplB protien of each organism. In order for MUSCLE to work, I have to put in an underscore on the fasta header.

The aligned file from Muscle in fasta format (rplB.afa). I use a python notebook (MakeTree.ipynb) to
make the phyloxml format. That file is used by graphlan to make the cool tree figure