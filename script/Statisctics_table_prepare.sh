sample="9,16,29,30,31,32,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,61,62,63,64,68,72,76,89,90,91,92,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,121,122,123,124,125"


###remove Archaea and mitochondia
cat result/table/otus_table.txt |grep -vE 'Archaea|mitochondria|Chlorophyta|Cyanobacteria'|sed '1d'|sed 's/#OTU ID/otu_id/g' > 1.data_for_statstics/otus_table_flt.txt

NF=`less 1.data_for_statstics/otus_table_flt.txt |awk -F"\t" '{print NF}'|head -1`
##remove sample from table
less 1.data_for_statstics/otus_table_flt.txt |csvtk transpose -t |csvtk sort -t  -k 1:n|csvtk grep -t -v -p "9,16,29,30,31,32,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,61,62,63,64,68,72,76,89,90,91,92,100,101,102,103,104,105,106,107,108,109,110,111,112,113,114,115,116,121,122,123,124,125"|csvtk transpose -t |csvtk cut -t -f 1,3-78,2 >1.data_for_statstics/otus_table_flt_sample.txt
#less 1.data_for_statstics/otus_table_flt.txt |csvtk transpose -t |csvtk sort -t  -k 1:n|csvtk grep -t -v -p $sample|csvtk transpose -t |csvtk cut -t -f 1,3-89,2 >1.data_for_statstics/otus_table_flt_sample.txt

cat 1.data_for_statstics/otus_table_flt_sample.txt|sed 's/otu_id/#OTU ID/g' >1.data_for_statstics/otus_table_group1.txt

biom convert -i 1.data_for_statstics/otus_table_group1.txt -o 1.data_for_statstics/otus_table_group1.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy

biom summarize-table -i 1.data_for_statstics/otus_table_group1.biom -o 1.data_for_statstics/group1_table_summary.txt
####multiple_rarefactions
multiple_rarefactions.py -i 1.data_for_statstics/otus_table_group1.biom -m 2300 -x 23000 -s 2300 -n 100 -o group1_rarefication --lineages_included &
###
biom convert -i group1_rarefication/rarefaction_23000_50.biom -o 1.data_for_statstics/group1_TIC.txt --to-tsv --header-key taxonomy

cat 1.data_for_statstics/group1_TIC.txt |sed '1d'|sed 's/#OTU ID/otu_id/g' > 1.data_for_statstics/group1_TIC_flt.txt

csvtk rename2 -t 1.data_for_statstics/group1_TIC_flt.txt -F -f '[0-9]*'  -p '(\d+)' -r 'X${1}'> 1.data_for_statstics/group1_TIC_flt_reheader.txt

biom convert -i 1.data_for_statstics/group1_TIC_flt_reheader.txt -o 1.data_for_statstics/group1_TIC_flt.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy

/home/galaxy/R/3.5.0/bin/Rscript script/ACM_table.R 1.data_for_statstics/group1_TIC_flt.txt 1.data_for_statstics/ACM_group1.txt

cat <(head -1 1.data_for_statstics/ACM_group1.txt|awk '{print "otu_id\t"$0}') <(sed '1d' 1.data_for_statstics/ACM_group1.txt) >1.data_for_statstics/tmp1.txt;mv 1.data_for_statstics/tmp1.txt 1.data_for_statstics/ACM_group1.txt

biom convert -i 1.data_for_statstics/ACM_group1.txt -o 1.data_for_statstics/group1_ACM.biom --to-hdf5 --table-type="OTU table" --process-obs-metadata taxonomy


sh script/chang_taxonomy_format.sh 1.data_for_statstics/ACM_group1.txt 1.data_for_statstics/ACM_group1_norm_tax_form.txt

#####beta_analysis
#---all_OTU-------
beta_diversity.py -i 1.data_for_statstics/otus_table_group1.biom -m weighted_unifrac -o 3.data_for_figure/w_uni_otu -t otus/rep_set.tre
#---TIC-------
beta_diversity.py -i 1.data_for_statstics/group1_TIC_flt.biom -m weighted_unifrac -o 3.data_for_figure/w_uni_tic -t otus/rep_set.tre
#---ACM----
beta_diversity.py -i 1.data_for_statstics/group1_ACM.biom -m weighted_unifrac -o 3.data_for_figure/w_uni_acm -t otus/rep_set.tre



###
cp 1.data_for_statstics/map.tsv 3.data_for_figure/
cp 1.data_for_statstics/ACM_group1_norm_tax_form.txt 3.data_for_figure/
/home/galaxy/R/3.5.0/bin/Rscript script/PL_R_script/1.Count_all_phyrum_family_RA.R /media/galaxy/My\ Passport/backup/olddata/03_ZHM/16s/shenglan/ otus_table_flt.txt ACM_group1_norm_tax_form.txt 5-6-7-8-9-20-25-30-35-36-37-38-39

##draw_fig_beta



/home/galaxy/R/3.5.0/bin/Rscript script/PL_R_script/3.Beta_diversity_phylo.R /media/galaxy/My\ Passport/backup/olddata/03_ZHM/16s/shenglan/ 3.data_for_figure/w_uni_tic/weighted_unifrac_group1_TIC_flt.txt Beta_diversity_TIC_w_unifrac.pdf 5-6-7-8-9-20-25-30-35-36-37-38-39 

/home/galaxy/R/3.5.0/bin/Rscript script/PL_R_script/3.Beta_diversity_phylo.R /media/galaxy/My\ Passport/backup/olddata/03_ZHM/16s/shenglan/ 3.data_for_figure/w_uni_acm/weighted_unifrac_group1_ACM.txt Beta_diversity_ACM_w_unifrac.pdf 5-6-7-8-9-20-25-30-35-36-37-38-39 
##draw_fig_sig_PCoA
/home/galaxy/R/3.5.0/bin/Rscript script/PL_R_script/5.significant_pcoa.R /media/galaxy/My\ Passport/backup/olddata/03_ZHM/16s/shenglan/ 5-6-7-8-9-20-25-30-35-36-37-38-39 2.stastics_analysis_result/shenglan_stat_otus_family_phyrum_result_with_FC_and_RA_table_group1.xlsx otu_rhizo_gp1_mutant_with_tax rhizo_sig.pdf Rhizo "gp1|5S|Col" 0.05,0.5,2
/home/galaxy/R/3.5.0/bin/Rscript script/PL_R_script/5.significant_pcoa.R /media/galaxy/My\ Passport/backup/olddata/03_ZHM/16s/shenglan/ 5-6-7-8-9-20-25-30-35-36-37-38-39 2.stastics_analysis_result/shenglan_stat_otus_family_phyrum_result_with_FC_and_RA_table_group1.xlsx otu_root_gp1_mutant_with_tax root_sig.pdf Root "gp1|35S|Col" 0.05,0.5,2

##draw_fig_all_otu_PCoA
/home/galaxy/R/3.5.0/bin/Rscript script/PL_R_script/7.PCOA_all_otu.R /media/galaxy/My\ Passport/backup/olddata/03_ZHM/16s/shenglan/ 5-6-7-8-9-20-25-30-35-36-37-38-39 rhizo.pdf Rhizo "gp1|35S|Col"
/home/galaxy/R/3.5.0/bin/Rscript script/PL_R_script/7.PCOA_all_otu.R /media/galaxy/My\ Passport/backup/olddata/03_ZHM/16s/shenglan/ 5-6-7-8-9-20-25-30-35-36-37-38-39 root.pdf Root "gp1|35S|Col"

##rarefication
/home/galaxy/R/3.5.0/bin/Rscript script/PL_R_script/4.FigS2.rarefaction_raw_otus.R /media/galaxy/My\ Passport/backup/olddata/03_ZHM/16s/shenglan/ 5-6-7-8-9-20-25-30-35-36-37-38-39 rare.pdf

#taxonomy

/home/galaxy/R/3.5.0/bin/Rscript script/PL_R_script/6.Taxonomy_overview.R /media/galaxy/My\ Passport/backup/olddata/03_ZHM/16s/shenglan/ 5-6-7-8-9-20-25-30-35-36-37-38-39 
