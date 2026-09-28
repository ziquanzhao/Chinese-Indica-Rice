#===============软件安装================
git clone https://github.com/Sh1ne111/KMERIA.git
cd KMERIA
conda env create -f kmeria_env.yaml
conda activate KMERIA

# Set the library path for htslib (adjust the path to your actual KMERIA directory)
export LD_LIBRARY_PATH=/your_path/KMERIA/lib:$LD_LIBRARY_PATH

# Change Permissions
chmod 755 /your_path/KMERIA/bin/* /your_path/KMERIA/external_tools/* /your_path/KMERIA/bimbamAsso/*

#Add KMERIA directories to your PATH environment variable
# KMERIA
export PATH="$PATH:/home/zzq/Software/KMERIA/bin/"
export PATH="$PATH:/home/zzq/Software/KMERIA/bimbamAsso/"
export PATH="$PATH:/home/zzq/Software/KMERIA/external_tools/"

conda install -c bioconda plink2 gemma gcta minimap2 spades
conda install bioconda::bowtie=1.3





#==========KMERIA获得kmer矩阵===============
#01_kmer_counts 为每个个体单独计数k-mer
kmeria count -k 31 -t 16 -m 5 -M 1000 -H 25-N2_k31_hist.txt -o ./01_kmer_counts/25-N2_k31.bin 25-N2_R1.fq.gz 25-N2_R2.fq.gz

#02_kmer_matrices 在种群层面创建k-mer计数矩阵
#sample_k31.list内容如下：
/mnt/g/Rice20260731/Kmer/01KmerCounts/25-N10_k31.bin
/mnt/g/Rice20260731/Kmer/01KmerCounts/25-N100_k31.bin
/mnt/g/Rice20260731/Kmer/01KmerCounts/25-N101_k31.bin
/mnt/g/Rice20260731/Kmer/01KmerCounts/25-N102_k31.bin
kmeria kctm -i ./XianSample_k31.bin.list -o ./XianSample_k31 -v --no-header -b 100000 -s 50000000 --buffer-size 5000000

#03_filtered_matrices k-mers过滤和校正
#生成样本测序深度和倍性文件
For Mixed ploidy species (XianSampleDepth.list):
25-N10	37.16719895	2
25-N100	37.44086842	2
25-N101	47.99431421	2
25-N102	36.47411447	2
25-N103	38.14483263	2
25-N104	32.52745974	2
25-N105	31.18304842	2

kmeria filter -i ./02_kmer_matrices -o ./03_filtered_matrices -d ./XianSampleDepth.list -t 16 -c 1000 -p 2 -s 0.6 --output-format compressed -v

# 04_kmer_to_BIMBAM 将k-mer计数矩阵转换为BIMBAM剂量格式
#二倍体
kmeria m2b --in ./03_filtered_matrices --out ./04_bimbam --input-format auto --threads 12 --bgzf-threads 6 --level 3 --buffer-size 256 --verbose --stats

#05_QK 基于 kmer 计算Q + K 计算
# kmer数目极多，进行PCA和亲缘关系分析时，可以通过随机采样来减少计算量，采样比例控制在0.1%~1%即可
# 使用kmeria sketch函数在第四步生成的每个bimabm.gz抽取一定数量的kmer，最后把结果追加到Sample_k31_sketch_n2000.bimbam中，这个就相当于LD剪枝之后的SNP.vcf.gz
for i in ../04KmerToBIMBAM/04_bimbam_616/*.bimbam.gz; do x=$(basename $i); echo "kmeria sketch -n 22000 $i > ${x%.bimbam.gz}_sketch_Totaln0.1.bimbam";done > 616Sample.sketch.sh
bash kmeria.sketch.sh
cat *sketch_Totaln0.1.bimbam > AllSample.k31.sketch_Totaln0.1.bimbam
#这里的-n 87000需要自己根据第四步获得的每个bimabm.gz文件中有多少个kmer，再根据0.1~1%比例去计算。
#计算Q+K
plink2 --vcf Sample_k31_sketch_n87000.vcf --make-bed --out Sample_k31_sketch_n87000 --allow-extra-chr
gcta64 --bfile Sample_k31_sketch_n87000 --make-grm --autosome-num 12 --out Sample_k31_sketch_n87000.GCTA
gcta64 --grm Sample_k31_sketch_n87000.GCTA --pca 100 --out Sample_k31_sketch_n87000.PCA
awk 'BEGIN{FS=" ";OFS="\t"} {print"1",$3,$4,$5,$6,$7}' Sample_k31_sketch_n87000.PCA.eigenvec > Sample_k31_sketch_n87000.PCA5.cov
gemma -bfile Sample_k31_sketch_n87000 -gk 2 -o Sample_k31_sketch_n87000.kinship.gk2 -outdir ./

# 06_KmerGWAS 基于 Kmer 进行 GWAS 分析
#使用gemma进行kmer GWAS分析
kmeria asso -i ../04KmerToBIMBAM/ -p PlantHeight.list -k ../05QK/Sample_k31_sketch_n87000.kinship.gk2.sXX.txt -c ../05QK/Sample_k31_sketch_n87000.PCA5.cov -o ./ -t 32 --bimbam-gzip --maf 0.01 --miss 0.1 --tool gemma -m lmm -n 1 --verbose
#-t 32 表示最多同时处理 32 个 BIMBAM 文件，也就是同时启动最多 32 个独立 GEMMA 进程，而不是让一个 GEMMA 进程使用 32 个线程
#PlantHeight.list列表如下，样本顺序注意和第五步SampleID.list一致，第一列是样本ID，第二列是表型. 参数-n 1，kmeria会自动+1,把第二列当作表型。这和原生gemma有差异。
25-N10	75.3
25-N100	82.6
25-N101	79.2
25-N102	81.7
25-N103	69.4
25-N104	68.9
25-N105	84.9

#07 GWAS 后将显著位点定位到基因组上
# 生成曼哈顿图时，我们不仅需要显著kmer，还需要一定量的非显著kmer，例如可以保留所有-log10(p)>=5的行，对于小于5的，可以抽取0.1%的kmer用于绘图
# 这里先仅仅保留-log10(p)>=5的行
for i in $(ls ../BLUE-Asso/*.bimbam.assoc.txt);do awk '$12 <= 1e-3' $i >> BLUE-PlantHeight.k31.gemma.3.list;done
awk '{printf ">BLUE-PlantHeight%d_%s\n%s\n", NR, $12, $2}' BLUE-GrainLength.k31.gemma.5.list > BLUE-GrainLength.k31.gemma.5.fa
#建立参考基因组索引
bowtie-build OsNIP.fa OsNIP
#比对
bowtie -f -v 1 --best --strata -m 1 -p 16 /mnt/g/Rice20260731/Kmer/06GWAS/OsNIP BLUE-GrainLength.k31.gemma.5.fa -S BLUE-GrainLength.k31.gemma.threshold5.v1.sam
#-f	输入FASTA	kmer一般是FASTA
#-v 1	1 mismatch	31 bp kmer要求完全匹配
#--best	搜索最佳比对	防止随机第一个位置
#--strata	只考虑最佳层	保证最佳匹配
#-m 1	只报告最多1个位置	唯一定位
#-p 16	多线程	加速
#从sam文件中解析出坐标
samtools view -F 4 BLUE-GrainLength.k31.gemma.threshold5.v0.sam | awk 'BEGIN{OFS="\t"; print "KmerID","Chr","Start","P-value"}{split($1,a,"_");print a[1],$3,$4,a[2]}' > BLUE-GrainLength.k31.gemma.threshold5.v0.Plot.list
#绘图
Rscript /mnt/g/Rice20260731/Kmer/06GWAS/KMERIA-GWAS.Plot.SignStarGene.R --input BLUE-PlantHeight.k31.gemma.3.bowtie.Final.Plot.list --star-gene StarGene_PlantHeight_candidates.xlsx --ld-annovation-gene 100000
