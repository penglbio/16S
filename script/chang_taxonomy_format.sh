col=`cat $1|head -1|awk -F"\t" '{print NF}'`
col1=`expr $col + 6`
less $1 |sed 's/; /\t/g'|sed 's/k__//g'|sed 's/p__//g'|sed 's/c__//g'|sed 's/o__//g'|sed 's/f__//g'|sed 's/g__//g'|sed 's/s__//g'|sed -E 's/\t+$//g'|sed 's/$/\t*\t*\t*\t*\t*\t*/g'|cut -f 1-$col1|sed 's/*/not_assigned/g'|sed -E 's/taxonomy\tnot_assigned\tnot_assigned\tnot_assigned\tnot_assigned\tnot_assigned\tnot_assigned/Domain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies/g'>$2
line=`cat $2|head -1 |awk '{print NF+1}'`
cut -f1-$line1 $2>${2}.txt
mv $2.txt $2
