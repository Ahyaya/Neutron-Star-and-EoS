#!/bin/bash
p_Loop=8
ResLine=$(cat Results.txt | wc -l)
p_Line=$((${ResLine}/${p_Loop}))
pf=1
while [ $pf -le ${p_Loop} ]
do
cat Results.txt | tail -n +$(((${pf}-1)*${p_Line}+1)) | head -n ${p_Line} > Results_p${pf}.txt
pf=$((${pf}+1))
done
echo "Part files deployed!"

