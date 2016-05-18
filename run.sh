#1->num test $2->flownum $3->rulenum $4->msizerate $5->ttlmax $6->interval $7runtime
for i in `seq 1 $1`
do
nohup ./build/autotest $i $2 $3 $4 $5 $6 $7 & 
done
