source activate qiime1
conda env export > doc/environment.yml

if [ ! -d fastqc ]; then
	mkdir -p cut_fastq fastqc/raw fastqc/clean clean stat assemble
fi
