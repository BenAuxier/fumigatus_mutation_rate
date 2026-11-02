REFERENCE=/lustre/BIF/nobackup/brigg002/novogene2024/Af293_combined.fna

##Here we will make the combined.vcf.gz that contains all four ancestors + progeny
#gatk CombineGVCFs -R $REFERENCE -O MR_AD.g.vcf.gz -V AD_ances_prog.list
#Joined genotyping - the combined variants will be spit out by this
#gatk IndexFeatureFile -I  MR_AD.g.vcf.gz
#gatk GenotypeGVCFs -R $REFERENCE -O MR_AD.vcf.gz -V MR_AD.g.vcf.gz

##Here we remove the isolates that are contamination (based on msh6 and TR)
#bcftools view -S ^contamination_AD.txt -Oz \
#        -o MR_AD_contamination_removed.vcf.gz MR_AD.vcf.gz

#echo "Start filtering DP"
#bcftools +setGT MR_AD_contamination_removed.vcf.gz \
#  -Oz -o MR_AD_DP.vcf.gz \
#  -- -tq -i 'FMT/DP<6' --new-gt .

#bcftools +setGT MR_AD_DP.vcf.gz \
#        -Oz -o MR_AD_DP_v02.vcf.gz \
#        -- -tq -i 'FMT/DP>0 & FMT/AD[:0]=0 & FMT/AD[:1]=0' --new-gt .

#bcftools +setGT MR_AD_DP_v02.vcf.gz \
#        -Oz -o MR_AD_DP_v03.vcf.gz \
#        -- -tq -i 'FMT/DP>8 & FMT/AD[:0]<5 & FMT/AD[:1]<5' --new-gt .

#bcftools index MR_AD_DP_v03.vcf.gz

gatk SplitVcfs -I MR_AD_DP_v03.vcf.gz --SNP_OUTPUT MR_AD_SNPs.vcf.gz --INDEL_OUTPUT MR_AD_INDELs.vcf.gz --STRICT false
echo "starting filtering SNPs"
#export PATH="/usr/bin"
bcftools filter -O u -s "QD" -m + -e "QD < 2.0" MR_AD_SNPs.vcf.gz | bcftools filter -O u -s "MQ" -m + -e "MQ < 40.0" |bcftools filter -O u -s "FS" -m + -e "FS > 60.0" | bcftools filter -O u -s "MQRank" -m + -e "MQRankSum < -12.5" | bcftools filter -O u -s "ReadPos" -m + -e "ReadPosRankSum < -8.0" | bcftools view -O z -f "PASS" > MR_AD_SNPs_filtered.vcf.gz
echo "finished filtering SNPs"
bcftools index MR_AD_SNPs_filtered.vcf.gz
bcftools view MR_AD_SNPs_filtered.vcf.gz | bcftools +fill-tags -Oz -o MR_AD_VAF_SNPs_filtered.vcf.gz -- -t FORMAT/VAF
bcftools index MR_AD_VAF_SNPs_filtered.vcf.gz
echo "VAF added to SNPs"
##Here we will adapt the VAF values
echo "Start filtering VAF in SNPs"
bcftools +setGT MR_AD_VAF_SNPs_filtered.vcf.gz \
  -Oz -o MR_AD_VAF_SNPs_adjusted.vcf.gz \
  -- -tq -i 'FMT/VAF>0.01 & FMT/VAF<0.99' --new-gt .

bcftools index MR_AD_VAF_SNPs_adjusted.vcf.gz
echo "Finished filtering VAF in SNPs"

echo "starting filtering INDELs"
bcftools filter -O u -s "QD" -m + -e "QD < 2.0" MR_AD_INDELs.vcf.gz | bcftools filter -O u -s "FS" -m + -e "FS > 200.0" | bcftools filter -O u -s "ReadPos" -m + -e "ReadPosRankSum < -20.0" | bcftools view -O z -f "PASS" > MR_AD_INDELs_filtered.vcf.gz
echo "finished filtering INDELs"
bcftools index MR_AD_INDELs_filtered.vcf.gz
bcftools view MR_AD_INDELs_filtered.vcf.gz | bcftools +fill-tags -Oz -o MR_AD_VAF_NDELs_filtered.vcf.gz -- -t FORMAT/VAF
bcftools index MR_AD_VAF_INDELs_filtered.vcf.gz
echo "VAF added to INDELs"
##Here we will adapt the VAF values

echo "Start filtering VAF in INDELs"
bcftools +setGT MR_AD_VAF_INDELs_filtered.vcf.gz \
  -Oz -o MR_AD_VAF_INDELs_adjusted.vcf.gz \
  -- -tq -i 'FMT/VAF>0.02 & FMT/VAF<0.98' --new-gt .

bcftools index MR_AD_VAF_INDELs_adjusted.vcf.gz
echo "Finished filtering INDELs in SNPs"
echo "If VAF shouldn't be included, GT has been set to NA (.)"
echo "Merge SNPs and INDELs"

bcftools concat -O z -a MR_AD_VAF_INDELs_adjusted.vcf.gz MR_AD_VAF_SNPs_adjusted.vcf.gz > MR_AD_VAF_adjusted.vcf.gz
gatk IndexFeatureFile -I MR_AD_VAF_adjusted.vcf.gz

echo "Start calculating % missing GT per variant"
bcftools +fill-tags -Oz -o MR_AD_missing_GT.vcf.gz MR_AD_VAF_adjusted.vcf.gz -- -t F_MISSING

echo "Start filtering variats with low % of GT"
bcftools view -i 'F_MISSING<=0.10' MR_AD_missing_GT.vcf.gz -Oz -o MR_AD_good_variants.vcf.gz
bcftools index MR_AD_good_variants.vcf.gz

echo "Start recalculating INFO field"
bcftools +fill-tags MR_AD_good_variants.vcf.gz -Oz -o tmp.vcf.gz -- -t AC,AN,AF
mv tmp.vcf.gz MR_AD_VAF_adjusted_v02.vcf.gz
bcftools index -f MR_AD_VAF_adjusted_v02.vcf.gz

echo "Start filtering on AC"
bcftools view -i "INFO/AC!=0 & INFO/AC!=INFO/AN" MR_AD_VAF_adjusted_v02.vcf.gz -Oz -o final_MR_AD_v02.vcf.gz
bcftools index -f final_MR_AD_v02.vcf.gz
echo "Done filtering on AC"

echo "Start subsampling for parents"
bcftools view -S AD_ances.list -Oz -o MR_AD_ances.vcf.gz final_MR_AD_v02.vcf.gz
bcftools index -f MR_AD_ances.vcf.gz
echo "Done subsampling parents"

echo "Start filtering on AC"
bcftools view -i "INFO/AC!=0 & INFO/AC!=4 & INFO/AN=4" MR_AD_ances.vcf.gz -Oz -o MR_AD_ances_v02.vcf.gz
bcftools index -f MR_AD_ances_v02.vcf.gz
echo "Done filtering on AC"

