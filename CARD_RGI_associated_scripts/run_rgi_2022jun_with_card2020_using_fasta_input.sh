#!/bin/bash

infile=$1 # Input file in fasta format, which is the required input format to CARD/RGI.
outdir=$2 # Output directory for the CARD/RGI output for one sample. Make sure each sample has its own output directory.

card_json_resource=/path/to/downloaded/card/json/file/card.json
rgi_singularity_container=./rgi_2022jun_singularity.sif

if [[ $card_json_resource == '/path/to/downloaded/card/json/file/card.json' ]]; then
  echo 'In script run_rgi_2022jun_with_card2020_using_fastq_input.sh, please set the variable $card_json_resource to where you have installed the CARD database file card.json.'
fi

currdir=$outdir
output_prefix="${currdir}"/output_rgi_2022jun_with_card2020_fasta_input
outfile="${output_prefix}"_output.txt

# The output directory will become the current directory when this script executes the command: cd $currdir
mkdir -p "${currdir}"
mkdir -p "${currdir}"/localDB

tmpdir="${outdir}"/tmp/output_rgi_2022jun_with_card2020_fasta_input
tmpdir_for_mplconfigdir="${outdir}"/tmp/MPLCONFIGDIR
mkdir -p $tmpdir
export MPLCONFIGDIR=$tmpdir_for_mplconfigdir

echo 'cd' $currdir
cd $currdir
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
echo '  --input_sequence' $infile '\'
echo '  --output_file' "${output_prefix}" '\'
echo '  --input_type contig --local --clean'

singularity run --cleanenv --no-home $rgi_singularity_container /rgi/rgi main \
  --input_sequence $infile \
  --output_file "${output_prefix}" \
  --input_type contig --local --clean
echo ''

echo 'Finished! outfile:' $outfile
echo ''

