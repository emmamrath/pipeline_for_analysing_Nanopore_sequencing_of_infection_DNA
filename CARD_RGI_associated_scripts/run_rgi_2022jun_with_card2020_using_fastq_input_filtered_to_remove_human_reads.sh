#!/bin/bash

infile=$1 # Input file in fastq.gz format. This script will reformat it as fasta format, for input to RGI.
in_centrifuge=$2 # Input file containing the Centrifuge _output_by_read output file, which identifies which sequencing reads are human DNA.
outdir=$3 # Output directory for the CARD/RGI output for one sample. Make sure each sample has its own output directory.

#card_json_resource=/path/to/downloaded/card/json/file/card.json
#singularity_directory=/path/to/downloaded/singularity/container
. ../config.sh

if [[ $card_json_resource == '/path/to/downloaded/card/json/file/card.json' ]]; then
  echo 'In script CARD_RGI_associated_scripts/run_rgi_2022jun_with_card2020_using_fastq_input_filtered_to_remove_human_reads.sh, please set the variable $card_json_resource to where you have installed the CARD database file card.json.'
fi
if [[ $singularity_directory == '/path/to/downloaded/singularity/containers' ]]; then
  echo 'In script CARD_RGI_associated_scripts/run_rgi_2022jun_with_card2020_using_fastq_input_filtered_to_remove_human_reads.sh, please set the variable $singularity_directory to where you downloaded the Singularity containers when you ran the command: wget https://github.com/emmamrath/pipeline_for_analysing_Nanopore_sequencing_of_infection_DNA/releases/download/Singularity_containers_v1/*.sif'
fi

rgi_singularity_container=$singularity_directory/rgi_2022jun_singularity.sif
htslib_singularity_container=$singularity_directory/htslib_singularity.sif

currdir=$outdir
output_prefix="${currdir}"/output_rgi_2022jun_with_card2020_filtered_to_remove_human_reads
outfile="${output_prefix}"_output.txt

# The output directory will become the current directory when this script executes the command: cd $currdir
mkdir -p "${currdir}"
mkdir -p "${currdir}"/localDB

tmpdir="${outdir}"/tmp/output_rgi_2022jun_with_card2020
tmpdir_for_mplconfigdir="${outdir}"/tmp/MPLCONFIGDIR
mkdir -p $tmpdir
export MPLCONFIGDIR=$tmpdir_for_mplconfigdir

in_fasta="${infile%.gz}"
tmp_fasta="${tmpdir}"/tmp_fasta.fasta
tmp_fasta_gz="${tmp_fasta}".gz
tmp_non_human="${tmpdir}"/non_human_reads.txt

# Using the Centrifuge file to identify which sequencing reads are human (taxid 9606).
echo 'cut -d$''\t' '-f1,3' $in_centrifuge '| grep -vP' '\t9606$' '| cut -d$''\t' '-f1 | sort -k1,1V | awk' '{print $0 " "}' '| uniq >' $tmp_non_human
cut -d$'\t' -f1,3 $in_centrifuge | grep -vP '\t9606$' | cut -d$'\t' -f1 | sort -k1,1V | uniq | awk '{print $0 " "}' > $tmp_non_human
echo ''

echo 'cd' $currdir
cd $currdir
echo ''

# Convert fastq input to fasta input format for rgi, at the same time filtering to remove human reads.
echo 'singularity run' $htslib_singularity_container 'bgzip -dc' $infile
singularity run $htslib_singularity_container bgzip -dc $infile
echo ''
echo 'awk' '{if ((NR % 4) == 1) {sub(/@/,">"); printf $0 "\t"}; if ((NR % 4) == 2) {print $0}}' $in_fasta '| \'
echo '  grep -f' $tmp_non_human '| sed -e' 's/\t/\n/g' '>' $tmp_fasta
awk '{if ((NR % 4) == 1) {sub(/@/,">"); printf $0 "\t"}; if ((NR % 4) == 2) {print $0}}' $in_fasta | \
  grep -f $tmp_non_human | sed -e 's/\t/\n/g' > $tmp_fasta
echo ''
echo 'singularity run' $htslib_singularity_container 'bgzip -f' $tmp_fasta
singularity run $htslib_singularity_container -f $tmp_fasta
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

rm -rf $tmp_fasta_gz $tmp_non_human
