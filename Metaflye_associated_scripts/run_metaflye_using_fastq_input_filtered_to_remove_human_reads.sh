#!/bin/bash

infile=$1 # Input file in fastq.gz format.
in_centrifuge=$2 # Input file containing the Centrifuge _output_by_read output file, which identifies which sequencing reads are human DNA.
outdir=$3 # Output directory for the Metaflye output. A directory below this directory will be created for this specific input, using the input file name.

metaflye_singularity_container=./metaflye_singularity.sif

htslib_singularity_container=$(dirname $metaflye_singularity_container)
htslib_singularity_container="${htslib_singularity_container}"/htslib_singularity.sif

indir=$(dirname $infile)
infile_basename=$(basename $infile)
infile_basename_prefix="${infile_basename%.gz}"
infile_basename_prefix="${infile_basename_prefix%.fastq}"
outdir_suffix=$infile_basename_prefix

tmpdir="${outdir}"/tmp/output_metaflye_using_fastq_input_filtered_to_remove_human_reads
mkdir -p $tmpdir

in_fastq="${infile%.gz}"
tmp_non_human="${tmpdir}"/"${sample}".non_human_reads.txt
tmp_fastq="${tmpdir}"/"${sample}".non_human_reads.fastq
tmp_fastq_gz="${tmp_fastq}".gz
tmp_fastq_gz_basename=$(basename $tmp_fastq_gz)
tmp_fastq_gz_dirname=$(dirname $tmp_fastq_gz)

# Using the Centrifuge file to identify which sequencing reads are human (taxid 9606).
echo 'cut -d$''\t' '-f1,3' $in_centrifuge '| grep -vP' '\t9606$' '| cut -d$''\t' '-f1 | sort -k1,1V | awk' '{print $0 " "}' '| uniq >' $tmp_non_human
cut -d$'\t' -f1,3 $in_centrifuge | grep -vP '\t9606$' | cut -d$'\t' -f1 | sort -k1,1V | uniq | awk '{print $0 " "}' > $tmp_non_human
echo ''

# Extract the fastq input that is not human
echo 'singularity run' $htslib_singularity_container 'bgzip -dc' $infile
singularity run $htslib_singularity_container bgzip -dc $infile
echo ''
echo 'awk '{'
echo '    if (((NR % 4) == 1) || ((NR % 4) == 2) || ((NR % 4) == 3)) {printf $0 "\t"}'
echo '    if ((NR % 4) == 0) {print $0}'
echo '  }' | \'
echo '  grep -f $tmp_non_human | sed -e' 's/\t/\n/g' '>' $tmp_fastq
awk '{
    if (((NR % 4) == 1) || ((NR % 4) == 2) || ((NR % 4) == 3)) {printf $0 "\t"}
    if ((NR % 4) == 0) {print $0}
  }' | \
  grep -f $tmp_non_human | sed -e 's/\t/\n/g' > $tmp_fastq
echo ''
echo 'singularity run' $htslib_singularity_container 'bgzip -f' $tmp_fastq
singularity run $htslib_singularity_container bgzip -f $tmp_fastq
echo ''

# Run Metaflye to construct potentially multiple contigs for multiple organisms from the input sequences.
echo 'singularity exec --bind' "${indir}"':/mnt1,'"${outdir}"':/mnt2 \'
echo '  '$metaflye_singularity_container '/home/Flye/bin/flye \'
echo '  --meta --nano-corr /mnt1/'"${infile_basename}" '--out-dir /mnt2/'"${outdir_suffix}"

singularity exec --bind "${tmp_fastq_gz_dirname}":/mnt1,"${outdir}":/mnt2 \
  $metaflye_singularity_container /home/Flye/bin/flye \
  --meta --nano-corr /mnt1/"${tmp_fastq_gz_basename}" --out-dir /mnt2/"${outdir_suffix}"
echo ''

echo 'Finished! outfile:' "${outdir}"/"${outdir_suffix}"/assembly.fasta
echo ''

