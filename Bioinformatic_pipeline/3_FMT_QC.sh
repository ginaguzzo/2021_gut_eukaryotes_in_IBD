#### Quality control pipline for FMT dataset ####

conda activate euks_in_ibd

mkdir -p fmt_data/raw_files/qc_reports

## Move raw files to directory
mv SRR1*.gz fmt_data/raw_files

## Generate quality control reports with fastQC and multiQC.
find fmt_data/raw_files/*.gz | parallel -j 10 -v 'fastqc {} -o fmt_data/raw_files/qc_reports'
multiqc fmt_data/raw_files/qc_reports -n multiqc_all_report -o fmt_data/raw_files/qc_reports/multiqc 
# The files are trimmed and processed.
# These have also been processed through KneadData, so will continue to marker gene analysis



###################
#### RAREFYING ####
###################

## Inspect read counts.
cat fmt_data/raw_files/qc_reports/multiqc/multiqc_all_report_data/multiqc_general_stats.txt | awk '{print $6}' | sort -n 

## Read counts: 1.7 to 32.1 million reads
## Rarefy to 1.7 million

## Sort fastq by headers before rarefying
mkdir -p fmt_data/rarefied_files/sorted_fastq/

find fmt_data/raw_files/*.gz | sed 's/.fastq.gz$//' | \
parallel -j 10 'zcat {}.fastq.gz | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > fmt_data/rarefied_files/sorted_fastq/{/}_sorted.fastq'

## Compress sorted fastqs
find fmt_data/rarefied_files/sorted_fastq/*.fastq | parallel -j 10 -v 'pigz --processes 2 --verbose {}' 

## Rarefy to 1.7 million reads
find fmt_data/rarefied_files/sorted_fastq/*.gz | sed 's/.fastq.gz$//' | \
parallel -j 10 'seqtk sample -s100 {}.fastq.gz 1700000 > fmt_data/rarefied_files/{/}_sub.fastq'

## Compress rarefied fastqs
find fmt_data/rarefied_files/*.fastq | parallel -j 10 -v 'pigz --processes 2 --verbose {}'

## Check rarefaction worked by counting reads
for file in $(find fmt_data/rarefied_files/*.fastq.gz) ; do \
echo $(zcat "$file" | wc -l)/4|bc ; \
done

## Check for errors in reads
for file in $(find fmt_data/rarefied_files/*.fastq.gz); do \
zcat -v "$file" | paste - - - - | awk -F"\t" '{ if (length($2) != length($4)) print $0 }' |  tr '\t' '\n' >> error_reads.txt; \
done
#If files appear in the above text file, redo rarefying step on those.

## After successfully rarefying, delete sorted reads.
rm -r fmt_data/rarefied_files/sorted_fastq
rm error_reads.txt





