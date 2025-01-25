#!/bin/bash

infile=$1 # Input fastq.gz file of sequenced DNA
outfile_prefix=$2 # Output files' prefix. The 2 output files will have _output_by_species.txt and _output_by_read.txt appended to the outfile_prefix respectively.

centrifuge_software=/path/to/where/one/has/installed/centrifuge
centrifuge_indexes=/path/to/where/one/has/installed/centrifuge/indexes

if [[ $centrifuge_software == '/path/to/where/one/has/installed/centrifuge' ]]; then
  echo 'In script run_centrifuge_on_one_fastq.sh, please set the variable $centrifuge_software to where you have installed the centrifuge software.'
fi
if [[ $centrifuge_indexes == '/path/to/where/one/has/installed/centrifuge/indexes' ]]; then
  echo 'In script run_centrifuge_on_one_fastq.sh, please set the variable $centrifuge_indexes to where you have installed the centrifuge reference data indexes for bacteria/viruses/etc.'
fi

export PATH="${centrifuge_software_path}"/install/bin:$PATH

echo 'centrifuge -x' "${centrifuge_indexes}"'/hpvc \'
echo '  -U' "${infile}" '\'
echo '  --report-file' "${outfile_prefix}"'_output_by_species.txt \'
echo '  -S' "${outfile_prefix}"'_output_by_read.txt'

centrifuge -x "${centrifuge_indexes}"/hpvc \
  -U "${infile}" \
  --report-file "${outfile_prefix}"_output_by_species.txt \
  -S "${outfile_prefix}"_output_by_read.txt
echo ''

echo 'Finished! outfile:' "${outfile_prefix}"_output_by_species.txt
echo '          outfile2:' "${outfile_prefix}"_output_by_read.txt
echo ''
