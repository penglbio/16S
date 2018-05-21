configfile: "config.yaml"

rule all:
	input:
		expand('barcode/{sample}_barcode.fq',sample=config['samples']),
		expand('cut_fastq/{sample}_R1.fastq.gz',sample=config['samples']),
		expand('cut_fastq/{sample}_R1.fastq.gz',sample=config['samples']),
		expand('clean/{sample}_R1_paired.fastq.gz', sample=config['samples']),
		expand('clean/{sample}_R2_paired.fastq.gz', sample=config['samples']),
		expand('fastqc/raw/{sample}_R1_fastqc.html', sample=config['samples']),
		expand('fastqc/raw/{sample}_R2_fastqc.html', sample=config['samples']),
		expand('fastqc/clean/{sample}_R1_paired_fastqc.html', sample=config['samples']),
		expand('fastqc/clean/{sample}_R2_paired_fastqc.html', sample=config['samples']), expand('stat/fastqc_stat.tsv'),
		expand('assemble/{sample}.extendedFrags.fastq',sample=config['samples']),
		expand('barcode/{sample}_barcode/barcodes.fastq',sample=config['samples']),
		expand('demuplix/{sample}/seqs.fna',sample=config['samples']),
		'demuplix/merge_seqs.fna',
		'split_seq',
		'chimeric/seqs_chimeric_filtered.fna',
		'otus',
		'result/table/otus_table.txt',
		'result/summary/taxa_summ',
		'result/table/rich_sparse_otu_table_summary.txt'

rule barcode_extract:
	input:
		r1='fastq/{sample}_R1.fastq.gz',
		r2='fastq/{sample}_R2.fastq.gz'
	output:
		barc='barcode/{sample}_barcode.fq',
		cuf1='cut_fastq/{sample}_R1.fastq.gz',
		cuf2='cut_fastq/{sample}_R2.fastq.gz'	
	shell:
		"paste <(seqkit subseq -r 5:8 {input.r1} |seqkit fx2tab ) <(seqkit subseq -r 5:8 {input.r2}|seqkit fx2tab)|csvtk cut -t -f 1,2,5,3,6|sed 's/\\t/****/'|sed 's/\\t//'|sed 's/\\t/****/'|sed 's/\\t//' |sed 's/\*\*\*\*/\\t/g'|sed -e 's/\:0\:[ATCG]*/:0:/'|awk -F'\\t' -vOFS='\\t' '{{print $1$2,$2,$3}}'|seqkit tab2fx >{output.barc} ; paste <(seqkit fx2tab {output.barc}|cut -f2) <(seqkit subseq -r 28:-1 {input.r1}|seqkit fx2tab)|sed -e 's/\:0\:[ATCG]*/:0:/'|awk -F'\\t' -vOFS='\\t' '{{print $2$1,$3,$4}}' |seqkit tab2fx |gzip -c > {output.cuf1} ; paste <(seqkit fx2tab {output.barc}|cut -f2) <(seqkit subseq -r 28:-1 {input.r2}|seqkit fx2tab)|sed -e 's/\:0\:[ATCG]*/:0:/'|awk -F'\\t' -vOFS='\\t' '{{print $2$1,$3,$4}}' |seqkit tab2fx |gzip -c > {output.cuf2}"

rule fastqc_raw_PE:
	input:
		'cut_fastq/{sample}_R1.fastq.gz',
		'cut_fastq/{sample}_R2.fastq.gz'
	output:
		'fastqc/raw/{sample}_R1_fastqc.html',
		'fastqc/raw/{sample}_R2_fastqc.html'
	shell:
		'fastqc -t 2 -o fastqc/raw {input}'

rule trimmomatic_PE:
	input:
		r1 = 'cut_fastq/{sample}_R1.fastq.gz',
		r2 = 'cut_fastq/{sample}_R2.fastq.gz'
	output:
		r1_paired = 'clean/{sample}_R1_paired.fastq.gz',
		r2_paired = 'clean/{sample}_R2_paired.fastq.gz',
		r1_unpaired = 'clean/{sample}_R1_unpaired.fastq.gz',
		r2_unpaired = 'clean/{sample}_R2_unpaired.fastq.gz'
	params:
		adapter = config['adapter']
	shell:
		'trimmomatic PE -threads 30 -phred33 {input.r1} {input.r2} {output.r1_paired} {output.r1_unpaired} {output.r2_paired} {output.r2_unpaired} ILLUMINACLIP:{params.adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:160'

rule fastqc_clean_PE:
	input:
		'clean/{sample}_R1_paired.fastq.gz',
		'clean/{sample}_R2_paired.fastq.gz'
	output:
		'fastqc/clean/{sample}_R1_paired_fastqc.html',
		'fastqc/clean/{sample}_R2_paired_fastqc.html'
	shell:
		'fastqc -t 2 -o fastqc/clean {input}'

rule fastqc_stat_PE:
	input:
		['fastqc/raw/{sample}_R1_fastqc.html'.format(sample=x) for x in config['samples']],
		['fastqc/raw/{sample}_R2_fastqc.html'.format(sample=x) for x in config['samples']],
		['fastqc/clean/{sample}_R1_paired_fastqc.html'.format(sample=x) for x in config['samples']],
		['fastqc/clean/{sample}_R2_paired_fastqc.html'.format(sample=x) for x in config['samples']]
	output:
		'stat/fastqc_stat.tsv'
	params:
		Rscript = config['Rscript_path']
	shell:
		'{params.Rscript} script/reads_stat_by_fastqcr.R'

rule sample_stat:
	input:
		r1='doc/{sample}_barcode.fa',
		r2='barcode/{sample}_barcode.fq',
		r3='doc/{sample}_barcode_table_seq'
	output:
		'stat/{sample}_sample_stat.tsv'
	shell:
		'script/sample_stat.sh {input.r1} {input.r2} {output} {input.r3}'

rule flash:
	input:
		r1='clean/{sample}_R1_paired.fastq.gz',
		r2='clean/{sample}_R2_paired.fastq.gz'
	output:
		out1='assemble/{sample}.extendedFrags.fastq',
		out2='assemble/{sample}'
	shell:
		'flash {input.r1} {input.r2} -o {output.out2} -t 30 2>&1 |tee assemble/flash.log ;touch {output.out2}'

rule extract_barcode:
	input:
		'assemble/{sample}.extendedFrags.fastq'
	output:
		o1='barcode/{sample}_barcode',
		o2='barcode/{sample}_barcode/barcodes.fastq'
	shell:
		"extract_barcodes.py -f {input} -o {output.o1} -c barcode_in_label --char_delineator '1:N:0:' -l 8"

rule split_library:
	input:
		r1='assemble/{sample}.extendedFrags.fastq',
		r2='barcode/{sample}_barcode/barcodes.fastq',
		map='doc/{sample}_map.tsv'	
	output:
		o1='demuplix/{sample}/seqs.fna',
		o2='demuplix/{sample}'
	shell:
		'split_libraries_fastq.py -o {output.o2} -i {input.r1} -b {input.r2} -m {input.map} --barcode_type 8 '	

rule merge_seq:
	input:
		['demuplix/{sample}/seqs.fna'.format(sample=x) for x in config['samples']]
	output:
		'demuplix/merge_seqs.fna'
	shell:
		'cat demuplix/*/seqs.fna > {output}'

rule split_file:
	input:
		'demuplix/merge_seqs.fna'
	output:
		'split_seq'
	shell:
		'seqkit split -p 128 -O split_seq {input}'	

rule remove_chimeric:
	input:
		r1='split_seq',
		r2='demuplix/merge_seqs.fna'
	params:
		ref=config['ref_gold_seq']	
	output:
		out1='chimeric',
		out2='chimeric/chimerics.txt',
		out3='chimeric/seqs_chimeric_filtered.fna'
	shell:
		"ls {input.r1}/merge_seqs.part_*|rush 'identify_chimeric_seqs.py -i {{}} -m usearch61 -o {output.out1}/{{%.}}_chimeric -r {params}' ; cat {output.out1}/*/chimeras.txt >{output.out2} ;filter_fasta.py -f {input.r2} -o {output.out3} -s {output.out2} -n"	

rule pick_otus:
	input:
		'chimeric/seqs_chimeric_filtered.fna'
	params:
		config['pick_otu_params']
	output:
		o1='otus',
		o2='otus/otu_table_mc2_w_tax_no_pynast_failures.biom'
	shell:
		'pick_open_reference_otus.py -i {input} -o {output.o1}  -p {params} --suppress_step4 --min_otu_size 2 -f '

rule otu_table:
	input:
		'otus/otu_table_mc2_w_tax_no_pynast_failures.biom'
	output:
		'result/table/otus_table.txt'
	shell:
		'biom convert -i {input} -o {output} --to-tsv --header-key taxonomy'

rule otu_summary:
	input:
		r1='otus/otu_table_mc2_w_tax_no_pynast_failures.biom',
		map='doc/map.tsv'	
	output:
		'result/summary/taxa_summ'	
	shell:
		'summarize_taxa_through_plots.py -i {input.r1} -o {output} -m {input.map} -s'

rule biom_summary:
	input:
		'otus/otu_table_mc2_w_tax_no_pynast_failures.biom'
	output:
		'result/table/rich_sparse_otu_table_summary.txt'
	shell:
		'biom summarize-table -i {input} -o {output}'
