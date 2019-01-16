## clean ##
#rm clean/*_unpaired.fastq.gz
#rm bam/*.bed*
#rm count/*.bam 
#rm igv.log 
#rm track/*.bedgraph 
cp -r fastqc result
cp -r stat result
cp -r arare result
cp -r bdiv_even result
pwd|sed 's/.*\///g'|rush 'tar cvzf {}.result.tar.gz result'
