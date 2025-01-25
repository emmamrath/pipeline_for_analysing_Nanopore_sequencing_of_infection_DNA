#!/bin/bash

infile_fastq=$1 # Input file of sequencing reads in fastq.gz format.
infile_draft_assembly=$2 # Input file of assembled contigs in fasta format. This input will be the assembly.fasta output file from Metaflye.
outdir=$3 # Output directory for the Medaka output. A directory below this directory will be created for this specific input sample, using the input file name.

medaka_singularity_container=./medaka_singularity.sif

indir_fastq=$(dirname $infile_fastq)
infile_basename=$(basename $infile_fastq)
infile_basename_prefix="${infile_basename%.gz}"
infile_basename_prefix="${infile_basename_prefix%.fastq}"
outdir_suffix=$infile_basename_prefix
outdir_for_sample="${outdir}"/"${outdir_suffix}"
indir_draft_assembly=$(dirname $infile_draft_assembly)
infile_draft_assembly_basename=$(basename $infile_draft_assembly)
outfile="${outdir_for_sample}"/medaka_consensus/consensus.fasta

NPROC=1
BASECALLS=/mnt1/"${infile_basename}"
DRAFT=/mnt2/"${infile_draft_assembly_basename}"
OUTDIR=/mnt_out/medaka_consensus

mkdir -p $outdir_for_sample

# Run Medaka to polish previously constructed contig assemblies that were assembled from Nanopore sequencing.
echo 'singularity exec --bind' "${indir_fastq}"':/mnt1,'"${indir_draft_assembly}"':/mnt2,'"${outdir_for_sample}"':/mnt_out \'
echo '  '$medaka_singularity_container '\'
echo '  env PATH=/medaka/venv/bin:'$PATH '/medaka/venv/bin/medaka_consensus -i' ${BASECALLS} '-d' ${DRAFT} '-o' ${OUTDIR} '-t' ${NPROC} '-m r941_min_high_g303'

singularity exec --bind "${indir_fastq}":/mnt1,"${indir_draft_assembly}":/mnt2,"${outdir_for_sample}":/mnt_out \
  $medaka_singularity_container \
  env PATH=/medaka/venv/bin:$PATH /medaka/venv/bin/medaka_consensus -i ${BASECALLS} -d ${DRAFT} -o ${OUTDIR} -t ${NPROC} -m r941_min_high_g303
echo ''

echo 'Finished! outfile:' $outfile
echo ''

