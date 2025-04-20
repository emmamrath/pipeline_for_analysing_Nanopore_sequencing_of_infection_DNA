#!/bin/bash
set -euo pipefail

sample=samples_identifier
in_fastq=/path/to/input/data/fastq/sample.fastq.gz
infile_centrifuge_species=/path/to/input/data/centrifuge/output/sample_output_by_species.txt
infile_centrifuge_reads=/path/to/input/data/centrifuge/output/sample_output_by_read.txt
indir_seqID_fasta=/path/to/input/data/fasta/downloaded/from/NCBI/by/etch_and_align_all_centrifuge_organisms_of_a_sample.script
outdir=/path/to/writable/output/directory/for/this/sample
tmpdir=/path/to/writable/temporary/directory/for/this/sample

#samtools_path=/path/to/samtools/software/directory
#minimap_path=/path/to/minimap/software/directory
#bwa_path=/path/to/bwa/software/directory
# Please note that this script uses Rscript and thus needs R and libraries to be installed.
. ../config.sh

if [[ $samtools_path == '/path/to/samtools/software/directory' ]]; then
  echo 'In script analysis_scripts/align_reads_to_centrifuge_organisms_having_100_or_more_reads.sh, please set the variable $samtools_path to where you have installed the samtools software.'
fi
if [[ $minimap_path == '/path/to/minimap/software/directory' ]]; then
  echo 'In script analysis_scripts/align_reads_to_centrifuge_organisms_having_100_or_more_reads.sh, please set the variable $minimap_path to where you have installed the minimap2 software.'
fi
if [[ $bwa_path == '/path/to/bwa/software/directory' ]]; then
  echo 'In script analysis_scripts/align_reads_to_centrifuge_organisms_having_100_or_more_reads.sh, please set the variable $bwa_path to where you have installed the bwa software.'
fi

outfile_genome="${outdir}"/"${sample}".genome.fasta
outfile_genome_gz="${outdir}"/"${sample}".gz
out_sam_unsorted="${outdir}"/"${sample}".out_sam_unsorted.sam
outfile_depth="${outdir}"/"${sample}".samtools_depth.txt
outfile_contigs="${outdir}"/"${sample}".contig_names.txt
outfile_depth_summary="${outdir}"/"${sample}".samtools_depth_summary.txt

# List the taxids from Centrifuge, so as to choose the taxids that Centrifuge matched to the most reads

tmpfile_taxids="${tmpdir}"/"${sample}".compare_centrifuge_num_reads_in_taxID_and_seqID.taxID_ge_100.txt
outfile_taxids="${outdir}"/"${sample}".taxids.sorted_by_numReadsInTaxID.txt

echo -e "name\ttaxID\ttaxRank\ttaxIDgenomeSize\tnumReadsInTaxID\tnumUniqueReadsInTaxID\tabundance\tseqID\tnumReadsInSeqID" > $outfile_sorted

# name                    taxID   taxRank genomeSize      numReads        numUniqueReads  abundance
# Shewanella              22      genus   5140018         1               0               0.0
# Myxococcus macrosporus  35      species 8973512         2               0               0.0
# Sorangium cellulosum    56      species 13907952        3               0               0.0

while read inline; do

  IFS=$'\t' read -r -a array <<< "$inline"
  name1="${array[0]}"
  taxID1="${array[1]}"
  taxRank1="${array[2]}"
  genomeSize1="${array[3]}"
  numReads1="${array[4]}"
  numUniqueReads1="${array[5]}"
  abundance1="${array[6]}"

  if [[ $name1 != "name" ]]; then
    if [[ $numReads1 -ge 100 ]]; then

      # readID          seqID           taxID   score   2ndBestScore    hitLength       queryLength     numMatches
      # ERR3077929.1    NZ_CP044230.1   2610896 2095    0               170             2901            1
      # ERR3077929.2    NZ_CP044230.1   2610896 8266    0               255             2097            1
      # ERR3077929.3    species         1850093 81      81              24              2073            3

      awk -v name="$name1" -v taxID="$taxID1" -v taxRank="$taxRank1" -v genomeSize="$genomeSize1" -v numReads="$numReads1" -v numUniqueReads="$numUniqueReads1" -v abundance="$abundance1" 'BEGIN {FS="\t";OFS="\t"} {
        if ($3 == taxID) {print name, taxID, taxRank, genomeSize, numReads, numUniqueReads, abundance, $2}}' $infile_centrifuge_reads | \
        sort | uniq -c | sed 's/  */ /g' | sed -e 's/^ //' | sed -e 's/ /\t/' | sort -k1,1Vr | awk 'BEGIN {FS="\t";OFS="\t"} {print $2,$3,$4,$5,$6,$7,$8,$9,$1}' >> $tmpfile_taxids

    fi
  fi
done < $infile_centrifuge_species

# sort by numReads (which is for taxID) so that most abundant taxID appears first
# then sort by numReadsInSeqID so that most abundant seqID appears first within all the seqID for each taxID
# for sorting, need to not have name as the first field because name can have spaces in it, which would cause the sort to pick the wrong field to sort on

echo 'grep -v' '^name' $tmpfile_taxids '| awk' 'BEGIN {FS="\t";OFS="\t"} {print $2,$3,$4,$5,$6,$7,$8,$9,$1}' '| sort -k4,4Vr -k8,8Vr | awk' 'BEGIN {FS="\t";OFS="\t"} {print $9,$1,$2,$3,$4,$5,$6,$7,$8}' '>>' $outfile_taxids
grep -v '^name' $tmpfile_taxids | awk 'BEGIN {FS="\t";OFS="\t"} {print $2,$3,$4,$5,$6,$7,$8,$9,$1}' | sort -k4,4Vr -k8,8Vr | awk 'BEGIN {FS="\t";OFS="\t"} {print $9,$1,$2,$3,$4,$5,$6,$7,$8}' >> $outfile_taxids
echo ''

# head $outfile_taxids from Centrifuge
# name    taxID   taxRank taxIDgenomeSize numReadsInTaxID numUniqueReadsInTaxID   abundance       seqID   numReadsInSeqID
# Enterobacteriaceae      543     family  10341043        212     4       0.0     family  212
# Citrobacter koseri      545     species 4735357 141     5       0.0     NZ_CP044097.1   24
# Citrobacter koseri      545     species 4735357 141     5       0.0     NZ_LR134204.1   24
# Citrobacter koseri      545     species 4735357 141     5       0.0     species 21
# Citrobacter koseri      545     species 4735357 141     5       0.0     NZ_CP022073.1   9
# Enterobacter    547     genus   11028548        117     2       0.0     genus   117
# Enterobacter cloacae    550     species 5356576 113     10      0.0     species 44
# Enterobacter cloacae    550     species 5356576 113     10      0.0     NZ_CP009852.1   21
# Enterobacter cloacae    550     species 5356576 113     10      0.0     NZ_CP009857.1   6
# Enterobacter cloacae    550     species 5356576 113     10      0.0     NZ_CP029720.1   5
# Escherichia coli        562     species 8437634 20326   6313    0.0     NZ_CP012379.1   1022
# Escherichia coli        562     species 8437634 20326   6313    0.0     NZ_CP031215.1   409
# Escherichia coli        562     species 8437634 20326   6313    0.0     NZ_CP033762.1   295

##### Create the reference fasta, consisting of the centrifuge organisms to which 100 or more sequencing reads were aligned
##### For each taxID (eg. E.coli is a taxID), include only the first species (eg. only NZ_CP012379.1, so don't include NZ_CP031215.1 and NZ_CP033762.1)
##### Stop including taxIDs after 100 taxIDs have been included. Hopefully all true infections and all plasmids will have been included

:>$outfile_genome
prev_taxID=''
num_taxID=0
max_num_taxID=100

echo 'Read in and choose the seqIDs:' $outfile_taxids
echo 'to create a genome for this sample:' $outfile_genome
echo ''

while read inline; do

  IFS=$'\t' read -r -a array <<< "$inline"
  name="${array[0]}"
  taxID="${array[1]}"
  taxRank="${array[2]}"
  taxIDgenomeSize="${array[3]}"
  numReadsInTaxID="${array[4]}"
  numUniqueReadsInTaxID="${array[5]}"
  abundance="${array[6]}"
  seqID="${array[7]}"
  numReadsInSeqID="${array[8]}"

  # If seqID from Centrifuge is actually one of the following instead of an actual seqID: class|family|genus|order|phylum|seqID|species|superkingdom
  # then there will not be fasta file for it. That is how such lines are ignored by this script.

  if [[ $taxID != "taxID" ]]; then
    if [[ $taxID != "9606" ]]; then
      if [[ $taxID != $prev_taxID ]]; then
        if [[ $num_taxID -le $max_num_taxID ]]; then
          num_taxID=$(( num_taxID + 1 ))
          if [[ $num_taxID -le $max_num_taxID ]]; then
            infile_seqID_fasta="${indir_seqID_fasta}"/"${seqID}".fasta.gz
            if [[ -f "$infile_seqID_fasta" ]]; then
              echo 'bgzip -dc' $infile_seqID_fasta '>>' $outfile_genome
              bgzip -dc $infile_seqID_fasta >> $outfile_genome
              prev_taxID=$taxID
              num_taxID=0
            fi
          fi
        else
          echo ''
          echo 'Fasta file does not exist and is needed to include in the reference genome:' $infile_seqID_fasta
          echo ''
        fi
      fi
    fi
  fi
done < $outfile_taxids
echo ''

echo 'outfile_genome' $outfile_genome
echo ''

echo 'bgzip -f' $outfile_genome
bgzip -f $outfile_genome
echo ''

echo $samtools_path'/samtools faidx' $outfile_genome_gz
$samtools_path/samtools faidx $outfile_genome_gz
echo ''

echo $bwa_path'/bwa index' $outfile_genome_gz
$bwa_path/bwa index $outfile_genome_gz
echo ''

##### Align this sample's sequencing reads to the reference fasta made specifically for this sample

# Increase I option to -I8g so as to avoid following error:
# [E::sam_parse1] no SQ lines present in the header
# --split-prefix is needed if the reference is large, and makes mappings more accurate than when not present and not all of the index is loaded into memory

echo $minimap_path'/minimap2 -ax map-ont' $outfile_genome_gz $in_fastq '-split-prefix >' $out_sam_unsorted
$minimap_path/minimap2 -ax map-ont $outfile_genome_gz $in_fastq -split-prefix > $out_sam_unsorted
echo ''

echo $samtools_path'/samtools sort' $out_sam_unsorted '-o' $outfile
$samtools_path/samtools sort $out_sam_unsorted -o $outfile
echo ''

echo $samtools_path'/samtools index' $outfile
$samtools_path/samtools index $outfile
echo ''

##### Do some coverage and depth calculations

echo $samtools_path'/samtools depth' $outfile '>' $outfile_depth
$samtools_path/samtools depth $outfile > $outfile_depth
echo ''

echo 'bgzip -dc' $outfile_genome_gz '| awk {output each contig seqID, name, length} >' $outfile_contigs
bgzip -dc $outfile_genome_gz | awk 'BEGIN {
  FS="\t";OFS="\t";id="";name="";seq_length=0
}
{
  char1=substr($1,1,1)
  if (char1==">") {
    if (id != "") { print id, name, seq_length }
    seq_length=0
    n=split($0,arr," ")
    id=substr(arr[1],2)
    name=arr[3]; for (i=4; i<=n; i++) {name=name" "arr[i]}
  } else {
    seq_length = seq_length + length
  }
} END {
  if (id!="") { print id, name, seq_length }
}' > $outfile_contigs
echo ''

echo 'Rscript calc_depth_by_samtools_for_one_sample.R' $outfile_depth $outfile_contigs $sample $outfile_depth_summary
Rscript calc_depth_by_samtools_for_one_sample.R $outfile_depth $outfile_contigs $sample $outfile_depth_summary
echo ''

echo 'Finished! outfile:' $outfile
echo '          outfile_genome:' $outfile_genome_gz
echo '          outfile_depth_summary:' $outfile_depth_summary
echo ''

rm -rf $out_sam_unsorted
#rm -rf  $outfile_depth $outfile_contigs
