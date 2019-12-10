


## Fetch data from FIMM cluster
## scRNAseq-data
me=$(whoami)
scrna_files=outs/filtered_gene_bc_matrices
local_folder=/Users/$me/Dropbox/gvhd_scrnaseq/data/scRNAseq
scrna_path=/fs/vault/fastq/fimm_ngs_mustjoki/5P_Pilot_FM2348/Count_H3GJ5DRXX

mkdir $local_folder/2013/
mkdir $local_folder/2017/

scp -r jhuuhtan@atlas.genome.helsinki.fi:$scrna_path/SI-GA-C1/$scrna_files $local_folder/SI-GA-C1/
scp -r jhuuhtan@atlas.genome.helsinki.fi:$scrna_path/SI-GA-C2/$scrna_files $local_folder/SI-GA-C2/

scp -r jhuuhtan@atlas.genome.helsinki.fi:$scrna_path/SI-GA-A3/$scrna_files $local_folder/SI-GA-A3/
scp -r jhuuhtan@atlas.genome.helsinki.fi:$scrna_path/SI-GA-A4/$scrna_files $local_folder/SI-GA-A4/
scp -r jhuuhtan@atlas.genome.helsinki.fi:$scrna_path/SI-GA-A5/$scrna_files $local_folder/SI-GA-A5/
scp -r jhuuhtan@atlas.genome.helsinki.fi:$scrna_path/SI-GA-A6/$scrna_files $local_folder/SI-GA-A6/
