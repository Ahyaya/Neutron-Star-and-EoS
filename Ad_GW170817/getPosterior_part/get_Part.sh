#!/bin/bash
p_Loop=5
ResLine=$(cat ReSamPara | wc -l)
p_Line=$((${ResLine}/${p_Loop}))
pf=1
while [ $pf -le ${p_Loop} ]
do
cat ReSamPara | tail -n +$(((${pf}-1)*${p_Line}+1)) | head -n ${p_Line} > ReSamPara_p${pf}
pf=$((${pf}+1))
done
echo "Part files deployed!"

