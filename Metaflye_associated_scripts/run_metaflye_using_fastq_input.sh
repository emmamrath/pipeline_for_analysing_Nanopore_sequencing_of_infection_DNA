#!/bin/bash

infile_fastq=$1 # Input file in fastq.gz format.
outdir=$2 # Output directory for the Metaflye output. A directory below this directory will be created for this specific input, using the input file name.

#singularity_directory=/path/to/downloaded/singularity/containers
. ../config.sh

if [[ $singularity_directory == '/path/to/downloaded/singularity/containers' ]]; then
  echo 'In script Metaflye_associated_scripts/Metaflye_associated_scripts/run_metaflye_using_fastq_input.sh, please set the variable $singularity_directory to where you downloaded the Singularity containers when you ran the command: wget https://github.com/emmamrath/pipeline_for_analysing_Nanopore_sequencing_of_infection_DNA/releases/download/Singularity_containers_v1/*.sif'
fi

metaflye_singularity_container=$singularity_directory/metaflye_singularity.sif

indir=$(dirname $infile_fastq)
infile_basename=$(basename $infile_fastq)
infile_basename_prefix="${infile_basename%.gz}"
infile_basename_prefix="${infile_basename_prefix%.fastq}"
outdir_suffix=$infile_basename_prefix
outfile="${outdir}"/"${outdir_suffix}"/assembly.fasta

# Run Metaflye to construct potentially multiple contigs for multiple organisms from the input sequences.
echo 'singularity exec --bind' "${indir}"':/mnt1,'"${outdir}"':/mnt2 \'
echo '  '$metaflye_singularity_container '/home/Flye/bin/flye \'
echo '  --meta --nano-corr /mnt1/'"${infile_basename}" '--out-dir /mnt2/'"${outdir_suffix}"

singularity exec --bind "${indir}":/mnt1,"${outdir}":/mnt2 \
  $metaflye_singularity_container /home/Flye/bin/flye \
  --meta --nano-corr /mnt1/"${infile_basename}" --out-dir /mnt2/"${outdir_suffix}"
echo ''

echo 'Finished! outfile:' $outfile
echo ''

