#### Quality control pipline for 500FG dataset ####

source activate euks_in_ibd

mkdir -p 500fg_data/raw_files/qc_reports

## Move raw files to directory
mv SRR5*.gz 500fg_data/raw_files

## Generate quality control reports with fastQC and multiQC.
find 500fg_data/raw_files/*.gz | parallel -j 10 -v 'fastqc {} -o 500fg_data/raw_files/qc_reports'
multiqc 500fg_data/raw_files/qc_reports -n multiqc_all_report -o 500fg_data/raw_files/qc_reports/multiqc 
#MultiQC report reveals that adapters are still present in the dataset.


###############
#### FASTP ####
###############

## Use fastp to trim adapters, polyG (-g) and polyX tails (-x), filter low complexity reads (-y), 
## reduce overrepresentation of reads (-p), and correct bases in overlapping regions (-c). 
mkdir -p 500fg_data/fastp_trimmed_data

find 500fg_data/raw_files/SRR5*_1.fastq.gz | sed -r 's/_1.fastq.gz$//' | \
parallel -j 5 -v 'fastp --thread 2 -g -x -y -p -V -c \
--in1 {}_1.fastq.gz \
--in2 {}_2.fastq.gz \
--out1 500fg_data/fastp_trimmed_data/{/}_1_paired_fastp.fastq.gz \
--out2 500fg_data/fastp_trimmed_data/{/}_2_paired_fastp.fastq.gz \
--json 500fg_data/fastp_trimmed_data/fastp.json \
--html 500fg_data/fastp_trimmed_data/fastp.html'

## Generate quality control reports to check fastp results.
mkdir 500fg_data/fastp_trimmed_data/qc_reports

find 500fg_data/fastp_trimmed_data/*.gz | parallel -j 10 -v 'fastqc {} -o 500fg_data/fastp_trimmed_data/qc_reports'
multiqc 500fg_data/fastp_trimmed_data/qc_reports -n multiqc_all_report -o 500fg_data/fastp_trimmed_data/qc_reports/multiqc 
#Fastp has successfully removed the adapters and improved the quality of the sequences.



###########################
#### FIX FASTQ HEADERS ####
###########################

## On the first iteration of this pipeline, the spaces in the fastq headers created difficulties for Kneaddata to match paired reads.
## Thus, the spaces will be replaced with colons before further processing.

## Unzip files.
gunzip -v 500fg_data/fastp_trimmed_data/*.fastq.gz

## Remove space in headers
for file in $(find 500fg_data/fastp_trimmed_data/*.fastq | sed -r 's/[.]fastq//' | uniq); do \
perl -pe 'if ($.%4==1){s/ /:/}' "$file".fastq > "$file"_tmp.fastq ; \
done

## Also replace period in headers with colon
for file in $(find 500fg_data/fastp_trimmed_data/*tmp.fastq | sed -r 's/tmp[.]fastq//' | uniq); do \
perl -pe 'if ($.%4==1){s/\./:/}' "$file"tmp.fastq > "$file"fixedheaders.fastq ; \
done

## Inspect updated headers and remove original and tmp files if successful.
cat 500fg_data/fastp_trimmed_data/SRR5127524_1_paired_fastp_fixedheaders.fastq | less #Inspect more files here as desired.
rm 500fg_data/fastp_trimmed_data/*tmp.fastq
rm 500fg_data/fastp_trimmed_data/*fastp.fastq

## Compress fastp output with pigz.
find 500fg_data/fastp_trimmed_data/*.fastq | parallel -j 4 -v 'pigz --best --processes 2 --verbose {}'



###################
#### KNEADDATA ####
###################

## Use Kneaddata to remove human contamination.
mkdir 500fg_data/kneaddata_output

find 500fg_data/fastp_trimmed_data/*_1_paired_fastp_fixedheaders.fastq.gz | sed 's/1_paired_fastp_fixedheaders.fastq.gz$//' | \
parallel -j 4 'kneaddata -v -t 2 \
--input {}1_paired_fastp_fixedheaders.fastq.gz --input {}2_paired_fastp_fixedheaders.fastq.gz \
--remove-intermediate-output --bypass-trim \
--reference-db tools/kneaddata/Homo_sapiens \
-o 500fg_data/kneaddata_output'

## Compress kneaddata output with pigz.
find 500fg_data/kneaddata_output/*.fastq | parallel -j 4 -v 'pigz --best --processes 2 --verbose {}'

## Sort kneaddata output.
mkdir 500fg_data/kneaddata_output/decontaminated
mkdir 500fg_data/kneaddata_output/contamination_and_unpaired
mkdir 500fg_data/kneaddata_output/log_files

## Calculate sizes of decontaminated reads that are unmatched.
ll -S 500fg_data/kneaddata_output/*unmatched_?.fastq.gz 

## The unmatched reads are quite small, most <1 mb so, while not ideal, we will discard them for downstream analyses.

mv 500fg_data/kneaddata_output/*paired_?.fastq.gz 500fg_data/kneaddata_output/decontaminated
mv 500fg_data/kneaddata_output/*unmatched_?.fastq.gz 500fg_data/kneaddata_output/contamination_and_unpaired
mv 500fg_data/kneaddata_output/*unmatched_?_contam.fastq.gz  500fg_data/kneaddata_output/contamination_and_unpaired
mv 500fg_data/kneaddata_output/*contam_?.fastq.gz 500fg_data/kneaddata_output/contamination_and_unpaired
mv 500fg_data/kneaddata_output/*.log 500fg_data/kneaddata_output/log_files


## Generate quality control reports to check kneaddata decontaminated results.
mkdir 500fg_data/kneaddata_output/decontaminated/qc_reports

find 500fg_data/kneaddata_output/decontaminated/*.gz | parallel -j 10 -v 'fastqc {} -o 500fg_data/kneaddata_output/decontaminated/qc_reports'
multiqc 500fg_data/kneaddata_output/decontaminated/qc_reports -n multiqc_all_report -o 500fg_data/kneaddata_output/decontaminated/qc_reports/multiqc 
#Sequence quality has been preserved even after kneaddata



###################
#### RAREFYING ####
###################

## To be able to compare abundances of taxa, the samples will need to be normalised by rarefying to the same
## number of reads for each sample.
## Rarefied data are only used with EukDetect comparing 1000IBD and 500FG hits.

## Inspect read counts.
cat 1000ibd_data/kneaddata_output/decontaminated/qc_reports/multiqc/multiqc_all_report_data/multiqc_general_stats.txt | awk '{print $6}' | sort -n 
cat 500fg_data/kneaddata_output/decontaminated/qc_reports/multiqc/multiqc_all_report_data/multiqc_general_stats.txt | awk '{print $6}' | sort -n 

## 1000IBD read counts: 2.7 to 25.7 million reads
## 500FG read counts: 2.3 to 33.2 million reads
## The lowest read count is 2.28 million in the 500FG dataset, so the samples in both datasets will be rarefied to that.

## Sort fastq by headers before rarefying
mkdir -p 500fg_data/kneaddata_output/decontaminated/rarefied/sorted_fastq/

find 500fg_data/kneaddata_output/decontaminated/*.gz | sed 's/.fastq.gz$//' | \
parallel -j 10 'zcat {}.fastq.gz | paste - - - - | sort -k1,1 -t " " | tr "\t" "\n" > 500fg_data/kneaddata_output/decontaminated/rarefied/sorted_fastq/{/}_sorted.fastq'

## Compress sorted fastqs
find 500fg_data/kneaddata_output/decontaminated/rarefied/sorted_fastq/*.fastq | parallel -j 10 -v 'pigz --processes 2 --verbose {}' 

## Rarefy to 2.28 million reads
find 500fg_data/kneaddata_output/decontaminated/rarefied/sorted_fastq/*.gz | sed 's/.fastq.gz$//' | \
parallel -j 10 'seqtk sample -s100 {}.fastq.gz 2280000 > 500fg_data/kneaddata_output/decontaminated/rarefied/{/}_sub.fastq'

## Compress rarefied fastqs
find 500fg_data/kneaddata_output/decontaminated/rarefied/*.fastq | parallel -j 10 -v 'pigz --processes 2 --verbose {}'

## Check rarefaction worked by counting reads
for file in $(find 500fg_data/kneaddata_output/decontaminated/rarefied/*.fastq.gz) ; do \
echo $(zcat "$file" | wc -l)/4|bc ; \
done

## Check for errors in reads
for file in $(find 500fg_data/kneaddata_output/decontaminated/rarefied/*.fastq.gz); do \
zcat -v "$file" | paste - - - - | awk -F"\t" '{ if (length($2) != length($4)) print $0 }' |  tr '\t' '\n' >> error_reads.txt; \
done
#If files appear in the above text file, redo rarefying step on those.

## After successfully rarefying, delete sorted reads.
rm -r 500fg_data/kneaddata_output/decontaminated/rarefied/sorted_fastq
rm error_reads.txt

## Rename files so that 1 and 2 are at the end
cd 500fg_data/kneaddata_output/decontaminated/rarefied/

for file in SRR*1_sorted_sub.fastq.gz; do mv "$file" "${file/1_sorted_sub/sorted_sub_1}"; done
for file in SRR*2_sorted_sub.fastq.gz; do mv "$file" "${file/2_sorted_sub/sorted_sub_2}"; done

cd ../../../