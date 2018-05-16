## clean ##
#rm clean/*_unpaired.fastq.gz
#rm bam/*.bed*
#rm count/*.bam 
#rm igv.log 
#rm track/*.bedgraph 
cp -r fastqc result
cp -r stat result
pwd|sed 's/.*\///g'|rush 'tar cvzf {}.result.tar.gz result'
