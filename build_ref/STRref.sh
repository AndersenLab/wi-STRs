
#https://github.com/HipSTR-Tool/HipSTR-references/blob/master/mouse/mouse_reference.md
#1, let's download this repository as it contains the scripts we need:

#git clone https://github.com/HipSTR-Tool/HipSTR-references.git
#cd HipSTR-references


# $1 :  folder name
# $2:   ref genome
# $3: gff


# Input parameters
OUTPUT_prefix_dir=$1
INPUT_genome=$2




##########################
## build STR reference ###
##########################

mkdir $OUTPUT_prefix_dir

cd $OUTPUT_prefix_dir

mkdir raw_fasta trf_results fixed_trf_results

cd raw_fasta

#2 split genome fasta by chr ,  install: https://github.com/mdshw5/pyfaidx/

 
faidx -x $INPUT_genome


cd ..


#3 We now use Tandem Repeats Finder to identify repeats on each chromosome:

for chrom in I II III IV V X MtDNA;
do
    echo raw_fasta/$chrom.fa trf_results 5
done | xargs -L 1 -P 30 ../HipSTR_ref/run_TRF.sh



#4 We then filter out repeats with a period longer than 6 and fix a few issues with incorrect TRF entries:

for chrom in I II III IV V X MtDNA;
do
    echo ../HipSTR_ref/fix_trf_output.py trf_results/$chrom.fa fixed_trf_results/$chrom.fa
done | xargs -L 1 -P 40 python



#5 We reformat the TRF entries and filter to only include repeats with a sufficiently high score:

files=""
for chrom in I II III IV V X MtDNA;
do
    files="$files,fixed_trf_results/$chrom.fa"
done
files=`echo $files | sed "s/,//"`
python ../HipSTR_ref/trf_parser.py $files > filtered_repeats.bed
bedtools sort -i filtered_repeats.bed > filtered_repeats.sorted.bed






#6 we merge overlapping STRs into single entries and filter repeats that fail merging:

python ../HipSTR_ref/analyze_overlaps.py filtered_repeats.sorted.bed pass.$OUTPUT_prefix_dir fail.$OUTPUT_prefix_dir



#7 We then remove any entries within 10bp of a failed merge region:

bedtools window -w 10 -a pass.$OUTPUT_prefix_dir -b fail.$OUTPUT_prefix_dir -v > pass.$OUTPUT_prefix_dir.r2


#8 To minimize the effects of nearby STRs on genotyping errors, we extract entries that aren't within 10bp of another entry or are within 10bp of one or more entries that all share the same period

bedtools merge -i pass.$OUTPUT_prefix_dir.r2 -c 4,6 -o collapse -d 10 | grep -v "," > pass.$OUTPUT_prefix_dir.r3

bedtools merge -i pass.$OUTPUT_prefix_dir.r2 -c 4,4,4,6 -o collapse,count_distinct,distinct,collapse -d 10 | grep "," | awk '$5 == 1' | awk -v OFS="\t" '{print $1, $2, $3, $6, $7}' | sed "s/,/\//g" >> pass.$OUTPUT_prefix_dir.r3




#9 Lastly, we construct the final reference for $1 

cat pass.$OUTPUT_prefix_dir.r3 | bedtools sort | awk -v OFS="\t" '{print $1, $2, $3, $4, ($3-$2+1)/$4, "STR_"NR, $5}' > $OUTPUT_prefix_dir.hipstr_reference.bed
 
cp $OUTPUT_prefix_dir.hipstr_reference.bed ../

# delete folder
cd ..
rm -rf $OUTPUT_prefix_dir






