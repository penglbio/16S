#!/bin/bash
####################################################################
#
#A (quite) simple submit script for a one or tow processor job
#
####################################################################
#
# SGE options
#
#Change to the current working directory upon starting of the job
#$ -cwd
#
# Specify the kind of shell script you use, for example, bash
#$ -S /bin/bash
#
# join the error and standard output streams
#$ -j y
#
#
# don't flood myself with e-mail
#$ -m e
#
# this is my e-mail address
##$ -M zihailing@live.cn
#
#where the format error go
#$ -e /psc/home/zihailing
#where the format output go
#$ -o /psc/home/zihailing
# notify me about pending SIG_STOP and SIG_KILL
#$ -notify
#
# Specify the array start ,end , step
#\$ -t 1-128:1 
# end of SGE stuff
#########################################################
# now execute my job:
mkdir test 0.data 1.fastqc 2.trim 3.merged 4.slout 5.nochimera 6.pick_otu 7.summary 8.result 9.anova_analysis

FILE=`ls test/*.R1.fastq.gz`
c=1
filelist[0]=head
filename[0]=head
for file in $FILE
	do
		files=`echo $file|sed -e 's/.R1.fastq.gz//g'`
		name=`echo $file |sed -E 's/test\/[0-9]*_(L[1-2]_[0-9]*)\..*/\1/g'`	
		#name=`echo $file|sed -e 's/\.\/test\/Sample.*\/HMZ.*_\([0-9][0-9]*\)_DNA.*_\(R[1-2]\)\(.*\)/S\1\3_\2/g'|sed 's/.fastq.gz//g'|sed 's/_001_R1//g'`
		filelist[$c]=$files
		filename[$c]=$name
		c=`expr $c + 1`
	done

echo ${#filelist[@]}
#step1 fastQC quality control
#zcat ${filelist[$SGE_TASK_ID]}.R1.fastq.gz > 0.data/${filename[$SGE_TASK_ID]}_R1.fastq
#zcat ${filelist[$SGE_TASK_ID]}.R2.fastq.gz > 0.data/${filename[$SGE_TASK_ID]}_R2.fastq
#fastqc 0.data/${filename[$SGE_TASK_ID]}_R1.fastq -o 1.fastqc
#fastqc 0.data/${filename[$SGE_TASK_ID]}_R2.fastq -o 1.fastqc
##step2 delete the lower quality reads
#java -jar /share/apps/prog/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 20 -phred33 -trimlog logfile 0.data/${filename[$SGE_TASK_ID]}_R1.fastq 0.data/${filename[$SGE_TASK_ID]}_R2.fastq 2.trim/${filename[$SGE_TASK_ID]}_R1_trim.fastq 2.trim/${filename[$SGE_TASK_ID]}_R1_unpair.fastq 2.trim/${filename[$SGE_TASK_ID]}_R2_trim.fastq  2.trim/${filename[$SGE_TASK_ID]}_R2_trim_unpair.fastq ILLUMINACLIP:/share/apps/prog/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:160
#step3 merges pair reads
#flash -r 251 -f 450 -s 45 -o 3.merged/${filename[$SGE_TASK_ID]} 2.trim/${filename[$SGE_TASK_ID]}_R1_trim.fastq 2.trim/${filename[$SGE_TASK_ID]}_R2_trim.fastq
#step4 get the barcode
#cat 3.merged/L1*extend* >4.slout/L1_forward_reads.fastq
#cat 3.merged/L2*extend* >4.slout/L2_forward_reads.fastq
#extract_barcodes.py -f 4.slout/forward_reads.fastq -o 4.slout/barcode -c barcode_in_label --char_delineator '1:N:0:' -l 16
####extract_barcodes.py -f 3.merged/${filename[$SGE_TASK_ID]}_11.extendedFrags.fastq -o 4.slout/barcode_${SGE_TASK_ID} -c barcode_in_label --char_delineator '1:N:0:' -l 16

#step5 delete the chimera
#split_libraries_fastq.py -o 4.slout/ -i 4.slout/L1_forward_reads.fastq,4.slout/L2_forward_reads.fastq -b 4.slout/barcode1/L1_barcode1.fastq,4.slout/barcode2/L2_barcodes.fastq -m L1_map1.txt,L2_map1.txt --barcode_type 16
#seqkit split 4.slout/seqs.fna -w 0 -O 5.nochimera -p 128
#identify_chimeric_seqs.py -i 5.nochimera/seqs.part$SGE_TASK_ID.fna -m usearch61 -o 5.nochimera/usearch_$SGE_TASK_ID -r gold.fa 
#cat 5.nochimera/usearch_*/chimeras.txt >5.nochimera/chimeras.txt
#filter_fasta.py -f 4.slout/seqs.fna -o 5.nochimera/seqs_chimeras_filtered.fna -s 5.nochimera/chimeras.txt -n
#step6 pick_otus
#pick_open_reference_otus.py -i 5.nochimera/seqs_chimeras_filtered.fna -o 6.pick_otu/otus -s 0.1 -p param.txt --suppress_step4 --min_otu_size 10
#cp 6.pick_otu/otus/otu_table_mc10_w_tax_no_pynast_failures.biom 7.summary/otus_table.biom
#biom convert -i 7.summary/otus_table.biom -o 7.summary/otus_table.txt --to-tsv --header-key taxonomy
#csvtk transpose -t -T 7.summary/otus_table.txt |csvtk sort -k1|csvtk transpose -t -T |less >7.summary/otus_table1.txt
#step7 plot must run in the GUI
#summarize_taxa_through_plots.py -o 7.summary/taxa_summary -i 7.summary/otus_table.biom -m map.txt -s
#step8 result and stastics
#cp -r 7.summary/* 8.result/
#rm 8.result/*.biom 
#rm log*.txt
#mkdir -p 7.summary/tmp
#seqkit stat 0.data/*_R1.fastq|sed -E 's/ +/\t/g'|awk -v OFS='\t' '{print $1,$4}'|sed 's/_R1.fastq\t/\t/g'|sed 's/0.data\///'>7.summary/tmp/1_1.txt
#seqkit stat 0.data/*_R2.fastq|sed -E 's/ +/\t/g'|awk -v OFS='\t' '{print $1,$4}'|sed 's/_R2.fastq\t/\t/g'|sed 's/0.data\///'>7.summary/tmp/1_2.txt
#seqkit stat 2.trim/*R1_trim.fastq|awk -v OFS='\t' '{print $1,$4}'|sed 's/.*\///g'|sed 's/_R1_trim.fastq//g'>7.summary/tmp/2_1.txt
#seqkit stat 2.trim/*R2_trim.fastq|awk -v OFS='\t' '{print $1,$4}'|sed 's/.*\///g'|sed 's/_R1_trim.fastq//g'>7.summary/tmp/2_2.txt
seqkit stat 3.merged/*.extendedFrags.fastq|awk -v OFS='\t' '{print $1,$4}' |sed 's/.*\///g'|sed 's/.extendedFrags.fastq//g'>7.summary/tmp/3.txt
cat <(echo -e "file\tseq_num") <(cat  5.nochimera/usearch_*/non_chim*|sed 's/\..*//g'|sort|uniq -c|awk -v OFS='\t' '{print $2,$1}')>7.summary/tmp/4.txt
cat <(echo -e "file\traw_R1\traw_R2\ttrim_R1\trim_R2\tmerged_reads\tdechimeras") <(csvtk join -t 7.summary/tmp/*|sed '1d')>8.result/data_stastics.txt



#step9 Anova analysis
#min_number=`biom summarize-table -i 6.pick_otu/otus/otu_table_mc10_w_tax_no_pynast_failures.biom |grep 'Min:'|sed 's/Min://g'|sed 's/\..*//g'`
#let time=min_number/5
#multiple_rarefactions.py -i 6.pick_otu/otus/otu_table_mc10_w_tax_no_pynast_failures.biom -m $time -x $min_number -s $time -n 100 -o 9.anova_analysis/rarefied_otu_tables/
#alpha_diversity.py -i 9.anova_analysis/rarefied_otu_tables/ -m chao1,observed_otus,PD_whole_tree -o 9.anova_analysis/adiv_chao1_observed_otus_pd/ -t 6.pick_otu/otus/rep_set.tre
#biom convert -i 9.anova_analysis/rarefied_otu_tables/rarefaction_${min_number}_32.biom -o 9.anova_analysis/TIC_otu.txt
#cat 9.anova_analysis/TIC_otu.txt|csvtk filter -f "2-31>20" --any -t >9.anova_analysis/ACMs_otu.txt
#python normlize.py
#cut -f 1-31 9.anova_analysis/normlize_primer_data.txt|csvtk transpose -t |sort -k1.2,1n|csvtk transpose -t>9.anova_analysis/ACMs_normalize_sorted.txt
#R CMD BATCH Anova_analysis.R














#end of job script
