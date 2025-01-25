#!/bin/bash

infile=$1 # Input file in fastq.gz format. This script will reformat it as fasta format, for input to RGI.
outdir=$2 # Output directory for the CARD/RGI output for one sample. Make sure each sample has its own output directory.

card_json_resource=/path/to/downloaded/card/json/file/card.json
rgi_singularity_container=./rgi_2023may_singularity.sif

dedupe_script=$(dirname $rgi_singularity_container)
dedupe_script="${dedupe_script}"/dedupe.sh

htslib_singularity_container=$(dirname $rgi_singularity_container)
htslib_singularity_container="${htslib_singularity_container}"/htslib_singularity.sif

if [[ $card_json_resource == '/path/to/downloaded/card/json/file/card.json' ]]; then
  echo 'In script run_rgi_2023may_nudge_cardl_using_fastq_input.sh, please set the variable $card_json_resource to where you have installed the CARD database file card.json.'
fi

currdir=$outdir
output_prefix="${currdir}"/output_rgi_2023may_main_cardl_fastq_input
outfile="${output_prefix}"_output.txt

# The output directory will become the current directory when this script executes the command: cd $currdir
mkdir -p "${currdir}"
mkdir -p "${currdir}"/localDB

tmpdir="${outdir}"/tmp/output_rgi_2023may_main_cardl_fastq_input
tmpdir_for_mplconfigdir="${outdir}"/tmp/MPLCONFIGDIR
mkdir -p $tmpdir
export MPLCONFIGDIR=$tmpdir_for_mplconfigdir

dedupe_fastq="${tmpdir}"/dedupe_fastq.fastq
dedupe_fastq_gz="${dedupe_fastq}".gz
tmp_fasta="$tmpdir}"/dedupe_fasta.fasta
tmp_fasta_gz="${tmp_fasta}".gz

echo 'cd' $currdir
cd $currdir
echo ''

# dedupe to remove fastq sequences when there are too many to input into rgi
echo $dedupe_script 'in='$in_fastq 'out='$dedupe_fastq 'overwrite=true'
$dedupe_script in=$in_fastq out=$dedupe_fastq overwrite=true
echo ''

# Convert fastq input to fasta input format for rgi
echo 'cat' $dedupe_fastq '| awk' '{if ((NR % 4) == 1) {sub(/@/,">"); print $0}; if ((NR % 4) == 2) {print $0}}' '>' $tmp_fasta
cat $dedupe_fastq | awk '{if ((NR % 4) == 1) {sub(/@/,">"); print $0}; if ((NR % 4) == 2) {print $0}}' > $tmp_fasta
echo ''

# Compress the fasta for input to RGI
echo 'singularity run' $htslib_singularity_container 'bgzip -f' $tmp_fasta
singularity run $htslib_singularity_container bgzip -f $tmp_fasta
echo ''

# RGI software requires a copy of the CARD data to exist in localDB in the current directory.
# Thus, create a softlink to the CARD data from localDB in the current directory.
echo 'rm -rf' "${currdir}"'/localDB/card.json'
rm -rf "${currdir}"/localDB/card.json
echo ''
echo 'ln -s' $card_json_resource "${currdir}"'/localDB/card.json'
ln -s $card_json_resource "${currdir}"/localDB/card.json
echo ''

# run rgi software using card2020 database, to identify antimicrobial mechanisms in fasta input DNA
echo 'singularity run --cleanenv --no-home' $rgi_singularity_container '/rgi/rgi main \'
echo '  --input_sequence' $tmp_fasta_gz '\'
echo '  --output_file' "${output_prefix}" '\'
echo '  --input_type contig --local --clean'

singularity run --cleanenv --no-home $rgi_singularity_container /rgi/rgi main \
  --input_sequence $tmp_fasta_gz \
  --output_file "${output_prefix}" \
  --input_type contig --local --clean
echo ''

echo 'Finished! outfile:' $outfile
echo ''

rm -rf $dedupe_fastq_gz $tmp_fasta_gz
