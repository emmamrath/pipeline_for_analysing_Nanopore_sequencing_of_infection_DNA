#!/bin/bash

infile_fasta=$1 # Input file of sequencing reads in fasta format. Recommend using Metaflye+Medaka output as input here.
outdir=$2 # Output directory for the Prokka output. A directory below this directory will be created for this specific input sample, using the input file name.

#singularity_directory=/path/to/downloaded/singularity/containers
. ../config.sh

if [[ $singularity_directory == '/path/to/downloaded/singularity/containers' ]]; then
  echo 'In script Prokka_associated_scripts/run_prokka_using_fasta_input.sh, please set the variable $singularity_directory to where you downloaded the Singularity containers when you ran the command: wget https://github.com/emmamrath/pipeline_for_analysing_Nanopore_sequencing_of_infection_DNA/releases/download/Singularity_containers_v1/*.sif'
fi

prokka_singularity_container=$singularity_directory/prokka_singularity.sif

indir_fasta=$(dirname $infile_fasta)
infile_basename=$(basename $infile_fasta)
infile_basename_prefix="${infile_basename%.gz}"
infile_basename_prefix="${infile_basename_prefix%.fasta}"
sample=$infile_basename_prefix
outdir_suffix=$infile_basename_prefix
outdir_for_sample="${outdir}"/"${outdir_suffix}"
outfile="${outdir_for_sample}"/"${sample}".faa
medaka_consensus/consensus.fasta

# Run Prokka
echo 'singularity exec' $prokka_singularity_container 'prokka' $infile_fasta '--outdir' $outdir_for_sample '--prefix' $sample '--force'
singularity exec $prokka_singularity_container prokka $infile_fasta --outdir $outdir_for_sample --prefix $sample --force
echo ''

echo 'Finished! outfile:' $outfile
echo ''

