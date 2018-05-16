cat <(seqkit grep -s -f $1 $2|seqkit seq -s|sort|uniq -c|sort -k1n|awk '{print $2"\t"$1}') <(seqkit grep -s -v -f doc/barcode.fa $2|seqkit seq -s|wc -l|awk '{print "nosample\t"$0}') >$3
Rscript script/sample_stat.R $3 $4 
