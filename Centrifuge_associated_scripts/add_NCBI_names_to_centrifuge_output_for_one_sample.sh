#!/bin/bash

infile=$1 # This is the _output_by_read.txt output from Centrifuge that contains sequencing read identifies and the NCBI nucleotide sequencing IDs that Centrifuge matched them to.

outfile="${infile%.txt}".plus_names.txt # This is the same as the input except that for every NCBI sequencing ID, there is also the NCBI sequence name so that we see what species it is.

tmpdir=$(dirname $outfile) # This script requires a temporary directory where files from NCBI can be downloaded, before being processed and deleted by this script.
sample=$(basename $outfile) # The temporary file is assigned a unique name in case this script is run for multiple sequences at the same time.
tmpdir="${tmpdir}"_"${sample}"
mkdir -p $tmpdir
tmp_list_ids="${tmpdir}"/tmp_list_ids.txt
tmp_list_ids_with_names="${tmpdir}"/tmp_list_ids_with_names.txt
tmp_fasta="${tmpdir}"/tmp_fasta.fasta

declare -a seqID_array
seqID_array=("class" "family" "genus" "kingdom" "no rank" "order" "phylum" "seqID" "species" "superkingdom" "unclassified")

echo 'Processing' $infile
echo ''

# get a list of the NCBI sequence IDs for each read
# but dont fetch human sequences of taxid 9606

echo 'cat' $infile '| awk' 'BEGIN {FS="\t"} {if ($3 != "9606") {print $0}}' '| cut -d$''\t' '-f2 | sort | uniq | grep -v seqID >' $tmp_list_ids
cat $infile | awk 'BEGIN {FS="\t"} {if ($3 != "9606") {print $0}}' | cut -d$'\t' -f2 | sort | uniq | grep -v seqID > $tmp_list_ids
echo ''

# for each NCBI sequence ID, get the organism name

:>$tmp_list_ids_with_names
while read seqID; do

  echo '...' $seqID
  organism_name="."
  if [[ ${seqID_array[*]} =~ $seqID ]]; then
    is_not_a_nucleotide_id=1
  else

    echo 'wget -O' $tmp_fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id=${seqID}"
    wget -O $tmp_fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id=${seqID}"
    head -n 1 $tmp_fasta
    organism_name=$(head -n 1 $tmp_fasta | sed 's/^>//')
  fi

  echo -e "${seqID}\t${organism_name}" >> $tmp_list_ids_with_names

done < $tmp_list_ids

echo 'Rscript add_NCBI_names_to_centrifuge_output_for_one_sample.R' $infile $tmp_list_ids_with_names $outfile
Rscript add_NCBI_names_to_centrifuge_output_for_one_sample.R $infile $tmp_list_ids_with_names $outfile
echo ''

echo 'Finished! outfile:' $outfile
echo ''

rm $tmp_list_ids $tmp_list_ids_with_names $tmp_fasta
