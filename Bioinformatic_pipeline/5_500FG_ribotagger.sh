#### RiboTagger for 500FG dataset ####


########################
#### RUN RIBOTAGGER ####
########################

mkdir -p 500fg_data/ribotagger/results

# Declare array and values
declare -a arr=("v4" "v5" "v6" "v7")

# Run RiboTagger for all 4 V regions #
for i in "${arr[@]}" 
do
    for file in $(ls 500fg_data/kneaddata_output/decontaminated/*.fastq.gz | sed -r 's/_[12][.]fastq[.]gz//' | uniq) 
    do
	tools/ribotagger-master/ribotagger/ribotagger.pl \
        -r "$i" \
        -no-bacteria -no-archaea \
        -in "${file}"_1.fastq.gz \
        "${file}"_2.fastq.gz \
        -out 500fg_data/ribotagger/results/"$(basename "$file")"."$i"
    done   
done


# RiboTagger BIOM files #

# Declare array and values
declare -a arr=("v4" "v5" "v6" "v7")

# RiboTagger biom reports for all 4 V regions     
for i in "${arr[@]}"
do 
    for file in $(ls 500fg_data/ribotagger/results/*"$i") 
    do 
        tools/ribotagger-master/ribotagger/biom.pl \
        -r "$i" \
        -in "${file}" \
        -out 500fg_data/ribotagger/results/"$(basename "$file")" \
        -taxonomy silva 
    done
done



#######################
#### FORMAT OUTPUT ####
#######################


# RiboTagger generating results table

mkdir 500fg_data/ribotagger/results/files_with_euks

cd 500fg_data/ribotagger/results

grep 'Eukaryota' ./*v?.anno >> files_with_euks.txt
sed -i -e 's/:.*//g' files_with_euks.txt # Remove end of line using wildcard option .* in sed
for file in $(cat files_with_euks.txt) ; do cp "$file" ./files_with_euks ; done
cd files_with_euks

# Declare array and values of V regions
declare -a arr=("v4" "v5" "v6" "v7")


for i in "${arr[@]}" 
do
    for file in $(find ./*"$i".anno)
    do 
        echo "$file" >> "$i"_hits.txt | sed -n "/Eukaryota/ s/$/\t$i/p" "$file" >> "$i"_hits.txt
    done
done

cat v4_hits.txt v5_hits.txt v6_hits.txt v7_hits.txt > all_ribotagger_v_region_hits.txt


# 500FG - format sample name per line
sed -i 's/_1_paired_fastp_fixedheaders_kneaddata_paired.v.*.anno//g' all_ribotagger_v_region_hits.txt
sed -i 's#./##g' all_ribotagger_v_region_hits.txt
awk -v OFS='\t' '{if($1 ~ "^SRR*") {tmp=$1} else {print tmp, $0}}' all_ribotagger_v_region_hits.txt > tmp.txt ; mv tmp.txt all_ribotagger_v_region_hits.txt


# Add header plus new column title for V region
cat *anno >> tmp.txt 
sed -n 1p tmp.txt > tmp2.txt ; mv tmp2.txt tmp.txt 
sed -i '1 s/$/\tv_region/' tmp.txt
sed -i '1 s/^/sampleid\t/' tmp.txt
cat all_ribotagger_v_region_hits.txt >> tmp.txt ; mv tmp.txt 500fg_ribotagger_all_v_region_hits.txt

# Copy results file to directory for R analysis
cp 500fg_ribotagger_all_v_region_hits.txt ../../../../../R_analysis/cohort_data_analysis/results
