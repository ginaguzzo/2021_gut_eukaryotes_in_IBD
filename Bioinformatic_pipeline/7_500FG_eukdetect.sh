#### EukDetect for 500FG dataset ####


#########################
#### UNRAREFIED DATA ####
#########################

mkdir -p 500fg_data/eukdetect/unrarefied_results

## First calculate the average read lengths of every sample, then the total average read length.
## This will be added to the eukdetect .yml file
cp calculate_*.sh 500fg_data/kneaddata_output/decontaminated/ 
cd 500fg_data/kneaddata_output/decontaminated/
bash calculate_avg_read_lengths_per_file.sh 
bash calculate_total_avg_read_lengths.sh
cd ../../../

## Edit configfile from default one.
#cp tools/EukDetect-master/default_configfile.yml tools/EukDetect-master/500fg_unrarefied_configfile.yml

## Total average is 80 bp, so this will be added to the config file.
## Keep paired end true and edit paths to directories.
## Remove sample1 and sample2 from sample list.
## Change _R1 and _R2 to _1 and _2.

## Add sample names and configure to be added to config file 
ls -1 500fg_data/kneaddata_output/decontaminated/*gz | xargs -n1 basename > 500fg_data/kneaddata_output/decontaminated/file_names_list.txt
sed -i 's/...........$//' 500fg_data/kneaddata_output/decontaminated/file_names_list.txt
sed -i 's/^/  /' 500fg_data/kneaddata_output/decontaminated/file_names_list.txt
sed -i 's/$/:/' 500fg_data/kneaddata_output/decontaminated/file_names_list.txt
cat 500fg_data/kneaddata_output/decontaminated/file_names_list.txt >> tools/EukDetect-master/500fg_unrarefied_configfile.yml

## Run eukdetect
source activate eukdetect

eukdetect --mode runall --configfile tools/EukDetect-master/500fg_unrarefied_configfile.yml --cores 10

mv snakemake_165216*log 500fg_data/eukdetect/


##################################
#### FORMAT OUTPUT - UNRAREFIED ##
##################################

cd 500fg_data/eukdetect/unrarefied_results

## Move files with hits to separate directory
mkdir files_with_euks
grep 'Taxid' ./*table.txt >> files_with_euks.txt
awk '{print $1}' files_with_euks.txt > tmp.txt ; mv tmp.txt files_with_euks.txt
sed -i 's/:Name//' files_with_euks.txt
sed 's/table/taxonomy/' files_with_euks.txt > tax_list.txt
cat tax_list.txt >> files_with_euks.txt ; rm tax_list.txt
for file in $(cat files_with_euks.txt) ; do cp "$file" ./files_with_euks ; done
cd files_with_euks

## Summarise samples with hits into one file
mkdir tables taxonomies
mv *table.txt tables
mv *taxonomy.txt taxonomies
cd tables
for file in $(ls ./*table.txt); do cat "$file"; done >> all_euk_results.txt
sed -n 1p all_euk_results.txt > tmp.txt
sed -i '/Name/d' all_euk_results.txt #remove line with Name in it
cat all_euk_results.txt >> tmp.txt ; mv tmp.txt all_euk_results.txt

## Format output file into usable table for R
awk -v OFS='\n' '{if($1=="Name"){print FILENAME} else {print $0}}' *table.txt > 500fg_eukdetect_all_hits_unrarefied.txt
sed -i 's/_1_paired_fastp_fixedheaders_kneaddata_paired_filtered_hits_table.txt//g' 500fg_eukdetect_all_hits_unrarefied.txt
# Remove EGA string at front of sample ID
awk -v OFS='\t' '{if($1 ~ "^SRR*") {tmp=$1} else {print tmp, $0}}' 500fg_eukdetect_all_hits_unrarefied.txt > tmp.txt ; mv tmp.txt 500fg_eukdetect_all_hits_unrarefied.txt

# Add header
for file in $(ls ./*table.txt); do cat "$file"; done >> tmp.txt
sed -ni 1p tmp.txt
cat 500fg_eukdetect_all_hits_unrarefied.txt >> tmp.txt ; mv tmp.txt 500fg_eukdetect_all_hits_unrarefied.txt
sed -i '1 s/^/sampleid\t/' 500fg_eukdetect_all_hits_unrarefied.txt

# Copy results file to directory for R analysis
cp 500fg_eukdetect_all_hits_unrarefied.txt ../../../../../../R_analysis/cohort_data_analysis/results


######################################################################################################################################################################


#######################
#### RAREFIED DATA ####
#######################

mkdir 500fg_data/eukdetect/rarefied_results

## First calculate the average read lengths of every sample, then the total average read length.
## This will be added to the eukdetect .yml file
cp calculate_*.sh 500fg_data/kneaddata_output/decontaminated/rarefied
cd 500fg_data/kneaddata_output/decontaminated/rarefied
bash calculate_avg_read_lengths_per_file.sh 
bash calculate_total_avg_read_lengths.sh
cd ../../../../

## Edit configfile from previous one.
#head -35 tools/EukDetect-master/500fg_unrarefied_configfile.yml > tools/EukDetect-master/500fg_rarefied_configfile.yml

## Total average is 88 bp, so this will be added to the config file.
## Keep paired end true and edit paths to directories.
## Remove sample1 and sample2 from sample list.
## Change _R1 and _R2 to _1 and _2.

## Add sample names and configure to be added to config file 
ls -1 500fg_data/kneaddata_output/decontaminated/rarefied/*gz | xargs -n1 basename > 500fg_data/kneaddata_output/decontaminated/rarefied/file_names_list.txt
sed -i 's/...........$//' 500fg_data/kneaddata_output/decontaminated/rarefied/file_names_list.txt
sed -i 's/^/  /' 500fg_data/kneaddata_output/decontaminated/rarefied/file_names_list.txt
sed -i 's/$/:/' 500fg_data/kneaddata_output/decontaminated/rarefied/file_names_list.txt
cat 500fg_data/kneaddata_output/decontaminated/rarefied/file_names_list.txt >> tools/EukDetect-master/500fg_rarefied_configfile.yml

## Run eukdetect
source activate eukdetect

eukdetect --mode runall --configfile tools/EukDetect-master/500fg_rarefied_configfile.yml --cores 10

mv snakemake_165216*log 500fg_data/eukdetect/


##################################
#### FORMAT OUTPUT - RAREFIED ####
##################################

## Formatting output into usable tables in R
cd 500fg_data/eukdetect/rarefied_results

## Move files with hits to separate directory
mkdir files_with_euks
grep 'Taxid' ./*table.txt >> files_with_euks.txt
awk '{print $1}' files_with_euks.txt > tmp.txt ; mv tmp.txt files_with_euks.txt
sed -i 's/:Name//' files_with_euks.txt
sed 's/table/taxonomy/' files_with_euks.txt > tax_list.txt
cat tax_list.txt >> files_with_euks.txt ; rm tax_list.txt
for file in $(cat files_with_euks.txt) ; do cp "$file" ./files_with_euks ; done
cd files_with_euks

## Summarise samples with hits into one file
mkdir tables taxonomies
mv *table.txt tables
mv *taxonomy.txt taxonomies
cd tables
for file in $(ls ./*table.txt); do cat "$file"; done >> all_euk_results.txt
sed -n 1p all_euk_results.txt > tmp.txt
sed -i '/Name/d' all_euk_results.txt #remove line with Name in it
cat all_euk_results.txt >> tmp.txt ; mv tmp.txt all_euk_results.txt

## Format output file into usable table for R
awk -v OFS='\n' '{if($1=="Name"){print FILENAME} else {print $0}}' *table.txt > 500fg_eukdetect_all_hits_rarefied.txt
sed -i 's/_1_paired_fastp_fixedheaders_kneaddata_paired_sorted_sub_filtered_hits_table.txt//g' 500fg_eukdetect_all_hits_rarefied.txt
# Remove EGA string at front of sample ID
awk -v OFS='\t' '{if($1 ~ "^SRR*") {tmp=$1} else {print tmp, $0}}' 500fg_eukdetect_all_hits_rarefied.txt > tmp.txt ; mv tmp.txt 500fg_eukdetect_all_hits_rarefied.txt

# Add header
for file in $(ls ./*table.txt); do cat "$file"; done >> tmp.txt
sed -ni 1p tmp.txt
cat 500fg_eukdetect_all_hits_rarefied.txt >> tmp.txt ; mv tmp.txt 500fg_eukdetect_all_hits_rarefied.txt
sed -i '1 s/^/sampleid\t/' 500fg_eukdetect_all_hits_rarefied.txt

# Copy results file to directory for R analysis
cp 500fg_eukdetect_all_hits_rarefied.txt ../../../../../../R_analysis/cohort_data_analysis/results

