#IBD分析
#生成一个默认1cM/Mb的遗传图谱，plink格式的，ChrID SNPid cM Pos
bcftools view -H 667XianSample.PASS.missing0.2.maf0.01.Biallelic.OnlyGT.OnlySNP.vcf.gz | awk 'BEGIN{OFS="\t"}{chrom=$1; pos=$2; id=$3; if(id=="."||id=="") id=chrom"_"pos;cm=pos/1e6; print chrom, id, cm, pos}' > 667XianSample.PASS.missing0.2.maf0.01.Biallelic.OnlyGT.OnlySNP.1cMperMb.map

# 运行hap-ibd
hap-ibd gt=667XianSample.PASS.missing0.2.maf0.01.Biallelic.OnlyGT.OnlySNP.vcf.gz map=667XianSample.PASS.missing0.2.maf0.01.Biallelic.OnlyGT.OnlySNP.1cMperMb.map out=667XianSample.PASS.missing0.2.maf0.01.Biallelic.OnlyGT.OnlySNP min-seed=0.5 min-extend=0.5 min-output=0.5 max-gap=1000 min-markers=180 min-mac=2 nthreads=24
#wget https://faculty.washington.edu/browning/ibd-ends.jar

#使用IBD.group.compara.py分析亚群间的IBD片段