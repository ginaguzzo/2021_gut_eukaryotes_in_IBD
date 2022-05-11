#### EukDetect for FMT dataset ####


#########################
#### UNRAREFIED DATA ####
#########################

mkdir -p fmt_data/eukdetect/unrarefied_results

## First calculate the average read lengths of every sample, then the total average read length.
## This will be added to the eukdetect .yml file
cp calculate_*.sh fmt_data/raw_files 
cd fmt_data/raw_files
bash calculate_avg_read_lengths_per_file.sh 
bash calculate_total_avg_read_lengths.sh
cd ../../

## Edit configfile from default one.
cp tools/EukDetect-master/default_configfile.yml tools/EukDetect-master/fmt_unrarefied_configfile.yml

## Total average is 130.42 bp, so 130 will be added to the config file.
## Change paired end to false and edit paths to directories.
## Remove sample1 and sample2 from sample list.


## Add sample names and configure to be added to config file 
ls -1 fmt_data/raw_files/*gz | xargs -n1 basename > fmt_data/raw_files/file_names_list.txt
sed -i 's/.........$//' fmt_data/raw_files/file_names_list.txt
sed -i 's/^/  /' fmt_data/raw_files/file_names_list.txt
sed -i 's/$/:/' fmt_data/raw_files/file_names_list.txt
cat fmt_data/raw_files/file_names_list.txt >> tools/EukDetect-master/fmt_unrarefied_configfile.yml

## Run eukdetect
source activate eukdetect

eukdetect --mode runall --configfile tools/EukDetect-master/fmt_unrarefied_configfile.yml --cores 10

mv snakemake_165216*log fmt_data/eukdetect/


##################################
#### FORMAT OUTPUT - UNRAREFIED ##
##################################

cd fmt_data/eukdetect/unrarefied_results

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
awk -v OFS='\n' '{if($1=="Name"){print FILENAME} else {print $0}}' *table.txt > fmt_eukdetect_all_hits_unrarefied.txt
sed -i 's/_filtered_hits_table.txt//g' fmt_eukdetect_all_hits_unrarefied.txt
# Remove EGA string at front of sample ID
awk -v OFS='\t' '{if($1 ~ "^SRR*") {tmp=$1} else {print tmp, $0}}' fmt_eukdetect_all_hits_unrarefied.txt > tmp.txt ; mv tmp.txt fmt_eukdetect_all_hits_unrarefied.txt

# Add header
for file in $(ls ./*table.txt); do cat "$file"; done >> tmp.txt
sed -ni 1p tmp.txt
cat fmt_eukdetect_all_hits_unrarefied.txt >> tmp.txt ; mv tmp.txt fmt_eukdetect_all_hits_unrarefied.txt
sed -i '1 s/^/sampleid\t/' fmt_eukdetect_all_hits_unrarefied.txt

# Copy results file to directory for R analysis
cp fmt_eukdetect_all_hits_unrarefied.txt R_analysis/fmt_data_analysis/results


###################################################################################################################################################################


#######################
#### RAREFIED DATA ####
#######################

mkdir fmt_data/eukdetect/rarefied_results

## First calculate the average read lengths of every sample, then the total average read length.
## This will be added to the eukdetect .yml file
cp calculate_*.sh fmt_data/rarefied_files 
cd fmt_data/rarefied_files
bash calculate_avg_read_lengths_per_file.sh 
bash calculate_total_avg_read_lengths.sh
cd ../../

## Edit configfile from default one.
head -35 tools/EukDetect-master/fmt_unrarefied_configfile.yml > tools/EukDetect-master/fmt_rarefied_configfile.yml

## Total average is 130.42 bp, so 130 will be added to the config file.
## Change paired end to false and edit paths to directories.
## Remove sample1 and sample2 from sample list.


## Add sample names and configure to be added to config file 
ls -1 fmt_data/rarefied_files/*gz | xargs -n1 basename > fmt_data/rarefied_files/file_names_list.txt
sed -i 's/.........$//' fmt_data/rarefied_files/file_names_list.txt
sed -i 's/^/  /' fmt_data/rarefied_files/file_names_list.txt
sed -i 's/$/:/' fmt_data/rarefied_files/file_names_list.txt
cat fmt_data/rarefied_files/file_names_list.txt >> tools/EukDetect-master/fmt_rarefied_configfile.yml

## Run eukdetect
source activate eukdetect

eukdetect --mode runall --configfile tools/EukDetect-master/fmt_rarefied_configfile.yml --cores 10

mv snakemake_165216*log fmt_data/eukdetect/


################################
#### FORMAT OUTPUT - RAREFIED ##
################################

cd fmt_data/eukdetect/rarefied_results

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
awk -v OFS='\n' '{if($1=="Name"){print FILENAME} else {print $0}}' *table.txt > fmt_eukdetect_all_hits_rarefied.txt
sed -i 's/_sorted_sub_filtered_hits_table.txt//g' fmt_eukdetect_all_hits_rarefied.txt
# Remove EGA string at front of sample ID
awk -v OFS='\t' '{if($1 ~ "^SRR*") {tmp=$1} else {print tmp, $0}}' fmt_eukdetect_all_hits_rarefied.txt > tmp.txt ; mv tmp.txt fmt_eukdetect_all_hits_rarefied.txt

# Add header
for file in $(ls ./*table.txt); do cat "$file"; done >> tmp.txt
sed -ni 1p tmp.txt
cat fmt_eukdetect_all_hits_rarefied.txt >> tmp.txt ; mv tmp.txt fmt_eukdetect_all_hits_rarefied.txt
sed -i '1 s/^/sampleid\t/' fmt_eukdetect_all_hits_rarefied.txt

# Copy results file to directory for R analysis
cp fmt_eukdetect_all_hits_rarefied.txt R_analysis/fmt_data_analysis/results

## There are only 13 samples with hits in the rarefied dataset, compared to 43 samples in the unrarefied dataset, so we will
## proceed with analysing only the unrarefied results in R.

