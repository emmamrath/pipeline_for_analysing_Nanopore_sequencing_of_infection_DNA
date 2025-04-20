#!/bin/bash
set -euo pipefail

sample=samples_identifier
infile_fastq=/path/to/input/data/fastq/sample.fastq.gz
infile_centrifuge_output=/path/to/input/data/centrifuge/output/sample_output_by_read.txt
outdir_centrifuge_output_plus_names=/path/to/writable/output/of/centrifuge/output/plus/organism/names
outdir_sample_genome=/path/to/writable/directory/where/a/reference/genome/for/this/sample/will/be/created # NCBI fasta files will be downloaded here too
outdir_sample_bam=/path/to/writable/directory/where/bam/alignment/file/will/be/created
tmpdir=/path/to/writable/temporary/directory/to/process/this/sample

#samtools_path=/path/to/samtools/software/directory
#minimap_path=/path/to/minimap/software/directory
# Please note that this script uses Rscript and thus needs R and libraries to be installed.
. ../config.sh

if [[ $samtools_path == '/path/to/samtools/software/directory' ]]; then
  echo 'In script analysis_scripts/fetch_and_align_all_centrifuge_organisms_of_a_sample.sh, please set the variable $samtools_path to where you have installed the samtools software.'
fi
if [[ $minimap_path == '/path/to/minimap/software/directory' ]]; then
  echo 'In script analysis_scripts/fetch_and_align_all_centrifuge_organisms_of_a_sample.sh, please set the variable $minimap_path to where you have installed the minimap2 software.'
fi

echo ''
echo '##################################################'
echo '# fetch_and_align_all_centrifuge_organisms_of_a_sample.sh'
echo '##################################################'
echo ''

echo '##################################################'
echo '# Processing sample:' $sample
echo '##################################################'
echo ''

echo ''
echo '##################################################'
echo '# add_names_to_output for sample' $sample
echo '##################################################'
echo ''

infile=$infile_centrifuge_output
indir=$(dirname $indir_centrifuge_output)
outdir=$outdir_centrifuge_output_plus_names

mkdir -p $tmpdir

tmp_list_ids="${tmpdir}"/tmp_list_ids.txt
tmp_list_ids_with_names="${tmpdir}"/tmp_list_ids_with_names.txt
tmp_fasta="${tmpdir}"/tmp_fasta.fasta
declare -a seqID_array
seqID_array=("class" "family" "genus" "kingdom" "no rank" "order" "phylum" "seqID" "species" "superkingdom" "unclassified")

echo 'Processing' $infile
echo ''
infile_basename=$(basename $infile)
outfile_basename="${infile_basename%.txt}".plus_names.txt
outfile="${outdir}"/"${outfile_basename}"
outfile_summary="${outdir}"/"${infile_basename%.txt}".plus_names_summary.txt

# get a list of the NCBI sequence IDs for each read

echo 'cat' $infile '| cut -d$''\t' '-f2 | sort | uniq | grep -v seqID >' $tmp_list_ids
cat $infile | cut -d$'\t' -f2 | sort | uniq | grep -v seqID > $tmp_list_ids
echo ''

# for each NCBI sequence ID, get the organism name

:>$tmp_list_ids_with_names
while read seqID; do

    echo '...' $seqID
    organism_name="."
    if [[ ${seqID_array[*]} =~ $seqID ]]; then
      is_not_a_nucleotide_id=1
    else

      # If this is a known human chromosome, don't need to fetch the sequence
      seqID_prefix=${seqID:0:7}
      if [[ $seqID_prefix == "NC_0000"  ]]; then

        IFS='.' read -r -a array <<< "$seqID_prefix"
        seqID_prefix_sans_version="${array[0]}"
        chrom=${seqID:7:2}
        if [ $chrom -lt 10 ]; then
          chrom=${chrom:2:1}
        fi
        if [[ $chrom == "23" ]]; then
          chrom="X"
        fi
        if [[ $chrom == "24" ]]; then
          chrom="Y"
        fi
        organism_name="${seqID} Homo sapiens chromosome ${chrom}"
        echo "This is a human chromosome, no need to fetch its name:" $organism_name

      else

        echo 'wget -O' $tmp_fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id=${seqID}" "> /dev/null 2>&1"
        wget -O $tmp_fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id=${seqID}" > /dev/null 2>&1
        #echo 'wget -O' $tmp_fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id=${seqID}" "-o /dev/null"
        #wget -O $tmp_fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id=${seqID}" -o /dev/null
        head -n 1 $tmp_fasta
        organism_name=$(head -n 1 $tmp_fasta | sed 's/^>//')

      fi
    fi

    echo -e "${seqID}\t${organism_name}" >> $tmp_list_ids_with_names

done < $tmp_list_ids

echo 'Rscript add_names_to_output.R' $infile $tmp_list_ids_with_names $outfile $outfile_summary
Rscript add_names_to_output.R $infile $tmp_list_ids_with_names $outfile $outfile_summary
echo ''

echo 'add_names_to_output outfile:' $outfile
echo 'add_names_to_output outfile_summary:' $outfile_summary
echo ''

echo ''
echo '##################################################'
echo '# top_list_ncbi_nucleotide_ids for' $sample
echo '##################################################'
echo ''

sample_genome="${sample}"_top
infile=$outfile_summary # outfile from the previous section: add_names_to_output
outdir=$outdir_sample_genome
outfile="${outdir}"/"${sample_genome}"_list_ncbi_nucleotide_ids.txt
tmp_infile="${tmpdir}"/tmp_"${sample}"_output_by_read.plus_names_summary.txt

# taxID         seqID           organism_name                                                   num_reads       num_reads_in_taxID
# 2081703       NZ_CP027242.1   NZ_CP027242.1 Peptostreptococcaceae bacterium oral taxon 929 strain W2294 chromosome, complete genome   45      45
# 525919        NC_013171.1     NC_013171.1 Anaerococcus prevotii DSM 20548, complete sequence  43      43
# 9606          NC_000009.12    NC_000009.12 Homo sapiens chromosome 9, GRCh38.p14 Primary Assembly     5       43
# 334413        NC_010376.1     NC_010376.1 Finegoldia magna ATCC 29328, complete sequence      41      41

prev_word1=""
prev_word2=""

mkdir -p $outdir
:>$outfile

# $infile
# taxID   seqID   organism_name   num_reads       num_reads_in_taxID
# 9606    NW_021160021.1  NW_021160021.1 Homo sapiens chromosome 17 genomic patch of type FIX, GRCh38.p14 PATCHES HG1320_PATCH    1       45850
# 9606    NW_021160025.1  NW_021160025.1 Homo sapiens chromosome 22 genomic patch of type FIX, GRCh38.p14 PATCHES HG494_PATCH     1       45850
# 158836  NZ_CP041054.1   NZ_CP041054.1 Enterobacter hormaechei strain C126 chromosome, complete genome   691     14489
# 158836  NZ_CP022532.1   NZ_CP022532.1 Enterobacter hormaechei strain MS7884A chromosome, complete genome        344     14489
# 158836  NZ_CP047570.1   NZ_CP047570.1 Enterobacter hormaechei strain F2 chromosome, complete genome     334     14489
# 158836  NZ_CP031565.1   NZ_CP031565.1 Enterobacter hormaechei strain 2013_1a chromosome, complete genome        330     14489

# $tmp_file
# taxID   seqID   organism_name   num_reads       num_reads_in_taxID
# 85698   NZ_LN890476.1   NZ_LN890476.1 Achromobacter xylosoxidans isolate R4 chromosome BN2905, complete sequence        1       1
# 1749    NZ_CP040635.1   NZ_CP040635.1 Acidipropionibacterium jensenii strain FAM 19038 chromosome, complete genome      1       3
# 1749    NZ_CP025570.1   NZ_CP025570.1 Acidipropionibacterium jensenii strain JS280 chromosome, complete genome  1       3
# 1749    NZ_LR134473.1   NZ_LR134473.1 Acidipropionibacterium jensenii strain NCTC13652 chromosome 1, complete sequence  1       3

echo 'Process' $infile
echo 'via' $tmp_infile
echo 'to produce' $outfile
echo ''

awk 'BEGIN {FS="\t";OFS="\t"} {
  n=split($3,a," ")
  new_name="."
  if (n > 1) {
    new_name=a[2]
    if (n > 2) {
      for (i=3; i<=n; i++) {
        new_name=new_name" "a[i]
      }
    }
  }
  print $5, $4, new_name, $0
}' $infile | sort -k3,3 -k4,4 -k1,1Vr -k2,2Vr | cut -d$'\t' -f4- > $tmp_infile

while read inline; do

  IFS=$'\t' read -r -a array <<< "$inline"
  taxID="${array[0]}"
  seqID="${array[1]}"
  organism_name="${array[2]}"
  num_reads="${array[3]}"
  num_reads_in_taxID="${array[4]}"

  ignore_this_record="0"
  if [[ $taxID == "taxID" ]]; then # ignore header
    ignore_this_record="1"
  fi
  if [[ $taxID == "9606" ]]; then # ignore human sequences
    ignore_this_record="1"
  fi
  if [[ $organism_name == "." ]]; then # ignore unclassified sequences
    ignore_this_record="1"
  fi
  if [[ $organism_name == "no rank ." ]]; then # ignore unclassified sequences
    ignore_this_record="1"
  fi

  if [[ $ignore_this_record == "0" ]]; then

      IFS=' ' read -r -a array2 <<< "$organism_name"
      # ignore the first word which is a repeat of the NCBI_id (array2 index 0)
      # Take the next two words (array2 indexes 1 and 2) as the short name of the species
      # for deciding in this loop whether we have seen this species already or not
      word1="${array2[1]}"
      word2="${array2[2]}"

      is_same_species="0"
      if [[ $prev_word1 == $word1 ]]; then
        if [[ $prev_word2 == $word2 ]]; then
          is_same_species="1"
        fi
      fi

      if [[ $is_same_species == "0" ]]; then

        # Rebuild the entire organism name as the original value in field 3
        # minus the first word which is a repeat of the NCBI id
        num_names="${#array2[@]}"
        org_name="${word1} ${word2}"
        for (( i=3; i<num_names; i++ )); do
          org_name="${org_name} ""${array2[i]}"
        done

        if [[ $org_name == "." ]]; then # ignore unclassified sequences
          ignore_this_record="1"
        fi
        if [[ $org_name == "" ]]; then # ignore unclassified sequences
          ignore_this_record="1"
        fi

        if [[ $ignore_this_record == "0" ]]; then

          echo -e "${seqID}\t${org_name}" >> $outfile

        fi

        prev_word1=$word1
        prev_word2=$word2

      fi

  fi

done < $tmp_infile

echo 'top_list_ncbi_nucleotide_ids outfile:' $outfile
echo ''

echo ''
echo '##################################################'
echo '# create_fasta_and_igv_genome_from_manual_fasta_for_one_sample for' $sample
echo '##################################################'
echo ''

sample_genome="${sample}"_top

infile_ids="${outdir_sample_genome}"/"${sample_genome}"_list_ncbi_nucleotide_ids.txt
outdir=$outdir_sample_genome
outfile_all_fasta="${outdir}"/"${sample_genome}".fasta
outfile_all_json="${outdir}"/"${sample_genome}".json
mkdir -p $outdir
tmpfile_seqid=$infile_ids

create_final_fasta=1
create_final_gff3=1
process_this_sample=1
if [[ -f "$outfile_all_fasta" ]]; then
  create_final_fasta=0
fi
if [[ -f "$outfile_all_json" ]]; then
  create_final_gff3=0
fi
if [[ $create_final_fasta -eq 0 ]]; then
  if [[ $create_final_gff3 -eq 0 ]]; then
    process_this_sample=0
  fi
fi

process_this_sample=1
if [[ $process_this_sample -eq 1 ]]; then

  # get NCBI sequences of organisms that centrifuge matches the reads to, they will be contigs for this sample

  ##### process each contig for this sample, by getting the fasta and genes in gff3

  while read nucleotide_id_from_centrifuge rest; do

    #IFS='.' read -r -a array <<< "$nucleotide_id_from_centrifuge"
    #nucleotide_id="${array[0]}"
    nucleotide_id=$nucleotide_id_from_centrifuge

    get_this_nucleotide=1

    if [[ $get_this_nucleotide -eq 1 ]]; then

      outfile_gff3="${outdir}"/"${nucleotide_id}".gff3
      outfile_fasta="${outdir}"/"${nucleotide_id}".fasta

      get_gff3=1
      get_fasta=1

      if [[ -f "$outfile_gff3" ]]; then
        file_size_gff3=$(stat -c %s $outfile_gff3)
        if [[ $file_size_gff3 -gt 0 ]]; then
          get_gff3=0
        fi
      fi
      if [[ -f "$outfile_fasta" ]]; then
        file_size_fasta=$(stat -c %s $outfile_fasta)
        if [[ $file_size_fasta -gt 0 ]]; then
          get_fasta=0
        fi
      fi

      if [[ $get_fasta -eq 1 ]]; then

        create_final_fasta=1

        #rm -rf $outfile_fasta
        wget -O $outfile_fasta "http://130.14.29.110/sviewer/viewer.cgi?db=nuccore&report=fasta&id=${nucleotide_id}"
        #wget -O $outfile_fasta "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=fasta&id=${nucleotide_id}"

        nucleotide_id_ncbi=$(head -n 1 $outfile_fasta | cut -d" " -f1 | sed 's/^>//')

        sed -i 's/^>/>'$nucleotide_id' /' $outfile_fasta
      fi

      if [[ $get_gff3 -eq 1 ]]; then

        create_final_gff3=1

        #rm -rf $outfile_gff3
        wget -O $outfile_gff3 "http://130.14.29.110/sviewer/viewer.cgi?db=nuccore&report=gff3&id=${nucleotide_id}"
        #wget -O $outfile_gff3 "https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi?db=nuccore&report=gff3&id=${nucleotide_id}"

        nucleotide_id_ncbi=$(head -n 1 $outfile_fasta | cut -d" " -f1 | sed 's/^>//')

        sed -i 's/'$nucleotide_id_ncbi'/'$nucleotide_id'/g' $outfile_gff3
      fi

    fi
  done < $tmpfile_seqid

  ##### create the final concatenated fasta and gff3 for this sample genome

  # create the beginning of the json file of gff3 gene definitions
  seen_first=0
  if [[ $create_final_gff3 -eq 1 ]]; then
    echo '{' > $outfile_all_json
    echo '  "id": "'$sample_genome'",' >> $outfile_all_json
    echo '  "name": "'$sample_genome'",' >> $outfile_all_json
    echo '  "fastaURL": "'$sample_genome'.fasta",' >> $outfile_all_json
    echo '  "indexURL": "'$sample_genome'.fasta.fai",' >> $outfile_all_json
    #echo '  "fastaURL": "'$outfile_all_fasta'",' >> $outfile_all_json
    #echo '  "indexURL": "'$outfile_all_fasta'.fai",' >> $outfile_all_json
    echo '  "tracks": [' >> $outfile_all_json
  fi

  # create the list of fasta files that will be used to create the multi-contig fasta for this sample
  list_fasta=""

  while read nucleotide_id_from_centrifuge rest; do

    #IFS='.' read -r -a array <<< "$nucleotide_id_from_centrifuge"
    #nucleotide_id="${array[0]}"
    nucleotide_id=$nucleotide_id_from_centrifuge
    echo 'Getting' $nucleotide_id

    get_this_nucleotide=1
    if [[ $get_this_nucleotide -eq 1 ]]; then

      outfile_gff3="${outdir}"/"${nucleotide_id}".gff3
      outfile_fasta="${outdir}"/"${nucleotide_id}".fasta

      # write the middle of the json file of gff3 gene definitions - one block per contig
      if [[ $create_final_gff3 -eq 1 ]]; then
        if [[ $seen_first -eq 1 ]]; then
          echo '    },' >> $outfile_all_json
        fi
        echo '    {' >> $outfile_all_json
        echo '      "name": "'$nucleotide_id'",' >> $outfile_all_json
        echo '      "format": "gff3",' >> $outfile_all_json
        echo '      "url": "'$nucleotide_id'.gff3"' >> $outfile_all_json
        #echo '      "url": "'$outfile_gff3'"' >> $outfile_all_json
        seen_first=1
      fi

      if [[ $create_final_fasta -eq 1 ]]; then
        if [[ $list_fasta == "" ]]; then
          list_fasta=$outfile_fasta
        else
          temp_list_fasta="${list_fasta} ${outfile_fasta}"
          list_fasta=$temp_list_fasta
        fi
      fi

    fi
  done < $tmpfile_seqid

  # write the end of the json file of gff3 gene definitions
  if [[ $create_final_gff3 -eq 1 ]]; then
    echo '    }' >> $outfile_all_json
    echo '  ]' >> $outfile_all_json
    echo '}' >> $outfile_all_json
  fi

  # cat a big fasta for this sample from the fasta of all its contigs
  if [[ $create_final_fasta -eq 1 ]]; then

    echo ''
    echo 'cat' $list_fasta '>' $outfile_all_fasta
    cat $list_fasta > $outfile_all_fasta
    echo ''

    echo $samtools_path'/samtools faidx' $outfile_all_fasta
    $samtools_path/samtools faidx $outfile_all_fasta
    echo ''

    echo 'rm' $list_fasta
    rm $list_fasta
    echo ''
  fi

fi

echo 'create_fasta_and_igv_genome_from_manual_fasta_for_one_sample outfile_all_fasta:' $outfile_all_fasta
echo 'create_fasta_and_igv_genome_from_manual_fasta_for_one_sample outfile_all_json:' $outfile_all_json
echo ''

echo ''
echo '##################################################'
echo '# align_sample_fastq_by_minimap2_to_manual_sample_genome_for_one_sample for' $sample
echo '##################################################'
echo ''

sample_fasta="${sample}"_top

NCPUS=1
ref_fasta="${outdir_sample_genome}"/"${sample_fasta}"/"${sample_fasta}".fasta
in_fastq=$infile_fastq
out_sam_unsorted="${outdir_sample_bam}"/"${sample_fasta}"_unsorted.sam
out_bam="${outdir_sample_bam}"/"${sample_fasta}".bam

# Didn't need this bwa index afterall
#echo bwa index' $ref_fasta
#bwa index $ref_fasta
#echo ''

# Increase I option to -I8g so as to avoid following error:
# [E::sam_parse1] no SQ lines present in the header

echo $minimap_path'/minimap2 -I 8G -ax map-ont' $ref_fasta $in_fastq '>' $out_sam_unsorted
$minimap_path/minimap2 -I 8G -ax map-ont $ref_fasta $in_fastq > $out_sam_unsorted
echo ''

echo $samtools_path'/samtools sort' $out_sam_unsorted '-o' $out_bam
$samtools_path/samtools sort $out_sam_unsorted -o $out_bam
echo ''

echo $samtools_path'/samtools index' $out_bam
$samtools_path/samtools index $out_bam
echo ''

echo 'align_sample_fastq_by_minimap2_to_manual_sample_genome_for_one_sample out_bam:' $out_bam
echo ''

echo ''
echo '##################################################'
echo '# calc_depth_by_samtools_for_one_sample for' $sample
echo '##################################################'
echo ''

bam_id="${sample}"_top

in_bam="${outdir_sample_bam}"/"${bam_id}".bam
in_fasta="${outdir_sample_genome}"/"${bam_id}"/"${bam_id}".fasta

outfile_depth="${in_bam%.bam}".samtools_depth.txt
outfile_contigs="${in_fasta%.fasta}".contig_names.txt
outfile_depth_summary="${in_bam%.bam}".samtools_depth_summary.txt

echo $samtools_path'/samtools depth' $in_bam '>' $outfile_depth
$samtools_path/samtools depth $in_bam > $outfile_depth
echo ''

# Don't use this code to get contigs info
#echo 'grep' '^>' $in_fasta '| cut -c2- | sed -e' 's/ /\t/' '| sed -e' 's/ /\t/' '| cut -d$''\t' '-f1,3- >' $outfile_contigs
#grep '^>' $in_fasta | cut -c2- | sed -e 's/ /\t/' | sed -e 's/ /\t/' | cut -d$'\t' -f1,3- > $outfile_contigs
#echo ''

# Instead use this code to get contigs info
echo 'awk' $in_fasta '>' $outfile_contigs
awk 'BEGIN {
  FS="\t";OFS="\t";id="";name="";seq_length=0
}
{
  char1=substr($1,1,1)
  if (char1==">") {
    if (id!="") { print id, name, seq_length }
    seq_length=0
    n=split($0,arr," ")
    id=substr(arr[1],2)
    name=arr[3]; for (i=4; i<=n; i++) {name=name" "arr[i]}
  } else {
    seq_length = seq_length + length
  }
} END {
  if (id!="") { print id, name, seq_length }
}' $in_fasta > $outfile_contigs
echo ''

echo 'Rscript calc_depth_by_samtools_for_one_sample.R' $outfile_depth $outfile_contigs $sample $outfile_depth_summary
Rscript calc_depth_by_samtools_for_one_sample.R $outfile_depth $outfile_contigs $sample $outfile_depth_summary
echo ''

echo 'num_unmapped=$('$samtools_path'/samtools view -f 4' $in_bam '| wc -l | cut -d' ' -f1 | awk ''{printf("%.0f\n", $0)}' ')'
num_unmapped=$($samtools_path/samtools view -f 4 $in_bam | wc -l | cut -d' ' -f1 | awk '{printf("%.0f\n", $0)}' )
echo ''

echo 'avg_length_unmapped=$('$samtools_path'/samtools view -f 4' $in_bam '| cut -d$''\t' '-f10 | awk' 'BEGIN {FS="\t";OFS="\t";sum=0} {sum=sum+length} END {avg_length=sum/NR; print avg_length}' '| cut -d' ' -f1 | awk ''{printf("%.0f\n", $0)}' ')'
avg_length_unmapped=$($samtools_path/samtools view -f 4 $in_bam | cut -d$'\t' -f10 | awk 'BEGIN {FS="\t";OFS="\t";sum=0} {sum=sum+length} END {avg_length=sum/NR; print avg_length}' | cut -d' ' -f1 | awk '{printf("%.0f\n", $0)}' )
echo ''

# The column headings of the output file from calc_depth_by_samtools_for_one_sample.R are:
# sample_id     ncbi_id contig_name     contig_length   avg_depth       num_nucleotides percent_nucleotides     avg_depth_ge10  num_nucleotides_depth_ge10  $
# num_unmapped will be the avg depth of unmapped sequences
# avg_length of unmapped sequences will be the num_nucleotides
# percent of contig unmapped is zero
n10=0
n100=0
n1000=0
a10=0
a100=0
a1000=0
if [[ $avg_length_unmapped -ge 10 ]]; then
  n10=$num_unmapped
  a10=$avg_length_unmapped
fi
if [[ $avg_length_unmapped -ge 100 ]]; then
  n100=$num_unmapped
  a100=$avg_length_unmapped
fi
if [[ $avg_length_unmapped -ge 1000 ]]; then
  n1000=$num_unmapped
  a1000=$avg_length_unmapped
fi

echo 'echo >>' $outfile_depth_summary
echo -e "${sample}\tunmapped\tunmapped\t0\t${num_unmapped}\t${avg_length_unmapped}\t0\t${n10}\t${a10}\t0\t${n100}\t${a100}\t0\t${n1000}\t${a1000}\t0" >> $outfile_depth_summary
echo ''

echo 'calc_depth_by_samtools_for_one_sample outfile_depth:' $outfile_depth
echo 'calc_depth_by_samtools_for_one_sample outfile_depth_summary:' $outfile_depth_summary
echo ''

echo ''
echo '##################################################'
echo '# collapse_depth_to_intervals for' $sample
echo '##################################################'
echo ''

infile=$outfile_depth
outfile="${outfile_depth%.txt}"_intervals.txt

# input is tab-delimited: contig   position   depth

echo 'awk -f' $infile '>' $outfile

awk 'BEGIN {
  FS="\t"; OFS="\t"
  curr_contig=""; curr_depth=0; curr_start=0; curr_end=0
}
{
  if (($1 == curr_contig) && ($3 == curr_depth)) {
    curr_end = $2
  } else {
    if (curr_depth > 0) {print curr_contig, curr_start, curr_end, curr_depth }
    curr_contig = $1; curr_start = $2; curr_end = $2; curr_depth = $3
  }
}
END {
  if (curr_depth > 0) {print curr_contig, curr_start, curr_end, curr_depth }
}' $infile > $outfile

echo 'collapse_depth_to_intervals outfile:' $outfile
echo ''

echo ''
echo '##################################################'
echo '# plot_the_coverage_for_all_contigs_in_a_sample for' $sample
echo '##################################################'
echo ''

infile_contigs=$outfile_contigs
infile_depth_intervals_file=$outfile
outfile_png="${infile_depth_intervals_file%.txt}".png

echo 'Rscript plot_the_coverage_for_all_contigs_in_a_sample.R' $sample $infile_contigs $infile_depth_intervals_file $outfile_png
Rscript plot_the_coverage_for_all_contigs_in_a_sample.R $sample $infile_contigs $infile_depth_intervals_file $outfile_png
echo ''

echo 'plot_the_coverage_for_all_contigs_in_a_sample outfile_png:' $outfile_png
echo ''

echo ''
echo '##################################################'
echo '# Finished!' $sample
echo '##################################################'
echo ''

