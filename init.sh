source activate qiime1
conda env export > doc/environment.yml

if [ ! -d fastqc ]; then
	mkdir -p doc cut_fastq fastqc/raw fastqc/clean clean stat assemble result/table result/summary 3.data_for_figure 1.data_for_statstics 2.stastics_analysis_result
fi
