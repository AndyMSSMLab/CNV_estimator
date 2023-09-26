VarIn=$1
PhenoIn=$2
CoVarIn=$3
Phenoname=$4
suffix=$5

regenie \
  --step 1 \
  --phenoCol ${Phenoname} \
  --pgen ${VarIn} \
  --bsize 100 \
  --lowmem --lowmem-prefix tmp_rg --loocv \
  --threads 18 \
  --bt \
  --spa \
  --phenoFile ${PhenoIn} \
  --covarFile ${CoVarIn} \
  --out Reg_S1_${suffix} 

regenie \
  --step 2 \
  --phenoCol ${Phenoname} \
  --pgen ${VarIn} \
  --bsize 100 \
  --bt \
  --spa \
  --phenoFile ${PhenoIn} \
  --covarFile ${CoVarIn} \
  --pThresh 0.05 \
  --pred Reg_S1_${suffix}_pred.list \
  --out Reg_S2_${foo}  
  

