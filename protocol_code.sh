########################################################################################################################
# Documents and scripts were written by: Jing Wu
# For manuscript: Jing Wu, et al. (2026). "Genomic origin, domestication, and gene expression control of fruit traits in Chinese bayberry with bHLH-IBH1 mediated anthocyanin biosynthesis" (under review). 
# email: 12307018@zju.edu.cn
# Jun Chen Lab. 
########################################################################################################################



########################################################################################################################
# Part 1 SNP identification in Myrica rubra
########################################################################################################################
## combine genotype gvcf
GenotypeGVCFs --output myrica_allSites.chr1.vcf.gz --include-non-variant-sites true --variant myrica_all.g.vcf.gz --intervals Chr1 --reference /mnt/d/cjevol/myrica/genome/Myrica_rubra_V2.fa
gatk MergeVcfs -I all_genotype_list.txt -O myrica_allSites.Genotype.vcf.gz
## call snp
gatk SelectVariants -V myrica_all.Genotype.vcf.gz -select-type SNP -O myrica_all.gatkRawSNP.vcf.gz
## filter
gatk VariantFiltration -R /data/wujing/Resequence/myrica/new_genome_mapping/ref/Myrica_rubra_V2.fa -V myrica_all.gatkRawSNP.vcf.gz -O myrica_all.gatkfilterSNP2.vcf.gz --filter-expression "QD < 2.5 || MQ < 40.0 || FS > 60.0 || SOR > 5.0 || QUAL < 30.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filter-name "Filter"
gatk SelectVariants -V myrica_all.gatkfilterSNP2.vcf.gz --exclude-filtered true -O myrica_all.filterd.vcf.gz
## Filter out low-quality samples
vcftools --gzvcf myrica_all.filterd.new.vcf.gz --recode --recode-INFO-all --remove remove25.txt > myrica_227.filterd.new.vcf.gz
## extract biallelic snp
vcftools --vcf myrica_227.filterd.new.vcf --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out myrica_227.2snp
## retain eight chromosomes
vcftools --vcf myrica_227.2snp.recode.vcf --chr Chr1 --chr Chr2 --chr Chr3 --chr Chr4 --chr Chr5 --chr Chr6 --chr Chr7 --chr Chr8 --recode --recode-INFO-all --out myrica_227.2snp.Chr8
## filter low DP
vcftools --vcf  myrica_227.2snp.Chr8.recode.vcf  --min-meanDP 5 --max-meanDP 50 --recode --out myrica_227.2snp.Chr8.DP
## fliter missing
vcftools --vcf  myrica_227.2snp.Chr8.DP.recode.vcf  --max-missing 0.95  --recode --out myrica_227.2snp.Chr8.DP.miss005
## fliter minor allele frequency
vcftools --vcf  myrica_227.2snp.Chr8.DP.miss005.recode.vcf --maf 0.05 --recode --out myrica_227.2snp.Chr8.DP.miss005.maf005
## filter LD
plink --chr-set 8 --bfile myrica_227.2snp.chr8.DP.miss005.maf005.recode --indep-pairwise 50 10 0.1 --out pruned
plink --chr-set 8 --bfile myrica_227.2snp.chr8.DP.miss005.maf005.recode --extract pruned.prune.in --make-bed --out myrica_227.snp.allfiltered.LD

########################################################################################################################
# Part 2 Neighbor-Joining (NJ) tree and ancestral genetic structures were inferred
########################################################################################################################
## Neighbor-Joining (NJ) tree
/data/wujing/software/VCF2Dis-1.53/bin/VCF2Dis  -InPut  /data/wujing/Resequence/myrica/new_genome_mapping/new0925/myrica_227.snp.allfiltered.LD.vcf  -OutPut p_227_snp_dis.mat
## ADMIXTURE
vcftools --vcf myrica_227.snp.allfiltered.LD.vcf --recode --recode-INFO-all --remove outgroup.txt > bayberry0.vcf
vcftools --vcf bayberry0.vcf --non-ref-ac-any 1 --recode --recode-INFO-all --out bayberry.vcf
vcftools --vcf bayberry.vcf.recode.vcf --plink --out bayberry
plink --chr-set 8 --file bayberry --make-bed --out bayberry
for K in 2 3 4 5 6 7 8 9 10 11 12 13 14 15; do /data/wujing/Resequence/myrica/new_genome_mapping/software/admixture/admixture_linux-1.3.0/admixture --cv ../bayberry.bed $K | tee log${K}.out; done
## Extract wild and cultivated individuals vcf
vcftools --vcf myrica_227.snp.allfilterd.LD.new.vcf --recode --recode-INFO-all --keep w1w2w3.txt > w1w2w3.snp.allfilterd.LD.new.vcf
vcftools --vcf myrica_277.snp.allfilterd.LD.new.vcf --recode --recode-INFO-all --keep cul.txt > cul.snp.allfilterd.LD.new.vcf
## PCA analysis
plink --bfile  /data/wujing/Resequence/myrica/new_genome_mapping/new0925/myrica_225.snp.allfiltered.LD.recode --allow-extra-chr --make-rel --pca 225 --out myrica_225.snp.allfiltered.LD.pca
plink --bfile  /data/wujing/Resequence/myrica/new_genome_mapping/new0925/w1w2w3.snp.allfilterd.LD.new.vcf --allow-extra-chr --make-rel --pca 131 --out w1w2w3.snp.allfilterd.LD.pca
plink --bfile  /data/wujing/Resequence/myrica/new_genome_mapping/new0925/cul.snp.allfilterd.LD.new.vcf --allow-extra-chr --make-rel --pca 131 --out cul.snp.allfilterd.LD.pca


########################################################################################################################
# Part 3  Genetic diversity parameters, LD decay and pairwise genetic correlation were calculated
########################################################################################################################
## genetic diversity parameters were calculated using pixy
bcftools concat --allow-overlaps Myrica_225_variant_filtered.vcf.gz  Myrica_225_invariant_filtered.vcf.gz -O z -o Myrica_225.allsites.miss90.vcf.gz
pixy --stats pi fst dxy --vcf Myrica_225.allsites.miss90.vcf.gz --populations popfile.txt --window_size 10000 --n_cores 64
## LD decay was calculated using PopLDdecay
/data/wujing/Resequence/myrica/new_genome_mapping/software/PopLDdecay-3.43/bin/PopLDdecay -InVCF  /data/wujing/Resequence/myrica/new_genome_mapping/xpehh/myrica_225.2snp.8.DP.miss005.maf005.vcf -OutStat  w1.stat.gz -SubPop w1.txt
/data/wujing/Resequence/myrica/new_genome_mapping/software/PopLDdecay-3.43/bin/PopLDdecay -InVCF  /data/wujing/Resequence/myrica/new_genome_mapping/xpehh/myrica_225.2snp.8.DP.miss005.maf005.vcf -OutStat  w2.stat.gz -SubPop w2.txt
/data/wujing/Resequence/myrica/new_genome_mapping/software/PopLDdecay-3.43/bin/PopLDdecay -InVCF  /data/wujing/Resequence/myrica/new_genome_mapping/xpehh/myrica_225.2snp.8.DP.miss005.maf005.vcf -OutStat  w3.stat.gz -SubPop w3.txt
/data/wujing/Resequence/myrica/new_genome_mapping/software/PopLDdecay-3.43/bin/PopLDdecay -InVCF  /data/wujing/Resequence/myrica/new_genome_mapping/xpehh/myrica_225.2snp.8.DP.miss005.maf005.vcf -OutStat  dk.stat.gz -SubPop dk.txt
/data/wujing/Resequence/myrica/new_genome_mapping/software/PopLDdecay-3.43/bin/PopLDdecay -InVCF  /data/wujing/Resequence/myrica/new_genome_mapping/xpehh/myrica_225.2snp.8.DP.miss005.maf005.vcf -OutStat  bq.stat.gz -SubPop bq.txt
/data/wujing/Resequence/myrica/new_genome_mapping/software/PopLDdecay-3.43/bin/PopLDdecay -InVCF  /data/wujing/Resequence/myrica/new_genome_mapping/xpehh/myrica_225.2snp.8.DP.miss005.maf005.vcf -OutStat  fh.stat.gz -SubPop fh.txt
/data/wujing/Resequence/myrica/new_genome_mapping/software/PopLDdecay-3.43/bin/PopLDdecay -InVCF  /data/wujing/Resequence/myrica/new_genome_mapping/xpehh/myrica_225.2snp.8.DP.miss005.maf005.vcf -OutStat  pl.stat.gz -SubPop pl.txt
perl /data/wujing/Resequence/myrica/new_genome_mapping/software/PopLDdecay-3.43/bin/Plot_MultiPop.pl -inList multi.list --output LD_fig -bin1 100 -bin2 1000
## kinship analysis
vcftools --vcf  myrica_225.2snp.8.DP.miss005.maf005.vcf --plink --out myrica_225.2snp.8.DP.miss005.maf005.recode
plink --file myrica_225.2snp.8.DP.miss005.maf005.recode --make-bed --out myrica_225.2snp.8.DP.miss005.maf005.recode
/data/wujing/software/king -b /data/wujing/Resequence/myrica/new_genome_mapping/xpehh/myrica_225.2snp.8.DP.miss005.maf005.recode.bed --kinship


########################################################################################################################
# Part 4    positive selection during domestication
########################################################################################################################
## XP-CLR analysis
for i in $(cat chr.txt);do xpclr --input /data/wujing/Resequence/myrica/new_genome_mapping/xpclr/myrica_227.2snp.DP.miss.maf.8.vcf --out /data/wujing/Resequence/myrica/new_genome_mapping/xpclr/"$i"_xpclr --samplesA c.txt --samplesB w.txt --chr "$i" --ld 0.95 --minsnps 20 --maxsnps 200 --size 20000 --step 2000 ;done
## XP-EHH analysis
## split eight chromosomes
touch 8chr.sh
#!bin/bash
for i in {1..8}
do
        vcftools --vcf myrica_225.2snp.8.DP.miss005.maf005.vcf --chr $i --recode --out $i
done
## replace 
touch findmiss.sh
#!/bin/bash
set -e
echo 'find missing'
echo $(date +%F%n%T)
for i in $(cat id.txt);do
         cat ./${i}.recode.vcf | perl -pe "s/\s\.:/\t.\/.:/g" | bgzip -c > ./${i}.findmissing.vcf.gz
         echo "${i} has been processed!"
done
echo 'Finished!'
echo $(date +%F%n%T)
## phasing
touch beagle.sh
#!/bin/bash
echo 'start phasing!' 
echo $(date +%F%n%T) 
for i in $(cat id.txt);do 
sample=${i}.findmissing.vcf.gz 
echo "Input ${i}.findmissing.gz file." 
echo $(date +%F%n%T) 
java -Xmx250g -jar beagle.22Jul22.46e.jar gt=${sample} out=${i}.phased nthreads=90 
echo "${i} has been phased!" 
echo $(date +%F%n%T) 
done 
echo 'Finished' 
echo $(date +%F%n%T)
## select wild individuals 
touch w.sh
#!/bin/bash
for i in $(cat id.txt );do
        echo "${i}.phased.vcf.gz is processing."
        echo $(date +%F%n%T)
        vcftools --gzvcf ${i}.phased.vcf.gz --keep ./w.txt --stdout --recode --recode-INFO-all | bgzip -c > ./w/${i}.w.phased.vcf.gz
        echo "${i}.w.phased.vcf.gz has been produced!"
        echo $(date +%F%n%T)
done
## select cultivated individuals 
touch c.sh
#!/bin/bash
for i in $(cat id.txt );do
        echo "${i}.phased.vcf.gz is processing."
        echo $(date +%F%n%T)
        vcftools --gzvcf ${i}.phased.vcf.gz --keep ./c.txt --stdout --recode --recode-INFO-all | bgzip -c > ./c/${i}.c.phased.vcf.gz
        echo "${i}.c.phased.vcf.gz has been produced!"
        echo $(date +%F%n%T)
done
## calculated
touch xpehh.sh
#!/bin/bash 
for i in {1..8};do 
/data/wujing/software/selscan-2.0.3/bin/linux/selscan --xpehh --vcf ${i}.c.phased.vcf.gz --vcf-ref ../wR/${i}.w.phased.vcf.gz --map ${i}.c.phased.map --pmap --threads 32 --out ./xpehh/${i}.cw.xpehh 
done
## normalization
touch norm.sh
#!/bin/bash 
for i in {1..8};do 
/data/wujing/software/selscan-2.0.3/bin/linux/norm --xpehh --files ${i}.cw.xpehh.xpehh.out --bp-win --winsize 100000 --log ${i}.cw.xpehh.norm.log --min-snps 20 
done
## combine
touch merge.sh
#!/bin/bash
for  ((i=1;i<=8;i++));do
        echo  cw.chr$i.xpehh.txt;done  |  xargs  -i  cat  {} >> cw.all.xpehh.txt
cat 0_header.txt cw.all.xpehh.txt > cw.final.xpehh.txt

########################################################################################################################
# Part 5    RNA data processing,mapping to reference genome using STAR,read count using featureCounts
########################################################################################################################
## mapping to reference genome
# step1 build index
STAR --runThreadN 6 --runMode genomeGenerate --genomeDir ./index_dir/ --genomeFastaFiles /data/wujing/Myrica/data/ref/Myrica_rubra_V2.fa --sjdbGTFfile /data/wujing/Myrica/data/ref/Myrica_rubra_V2.gtf --sjdbOverhang 149
# step2 run STAR
touch star.sh
#!bin/bash
for i in $(cat G1.txt);do STAR --twopassMode Basic --quantMode TranscriptomeSAM GeneCounts --runThreadN 20 --genomeDir /data/wujing/Myrica/mapping/index_dir --alignIntronMin 20 --alignIntronMax 50000 --outSAMtype BAM SortedByCoordinate --sjdbOverhang 149 --outSAMattrRGline ID:"$i" SM:"$i" PL:ILLUMINA --outFilterMismatchNmax 2 --outSJfilterReads Unique --outSAMmultNmax 1 --outFileNamePrefix /data/wujing/Myrica/mapping/$i/$i --outSAMmapqUnique 60 --readFilesCommand gunzip -c --readFilesIn /data/wujing/Myrica/data/fruit-A/"$i"_1.clean.fq.gz /data/wujing/Myrica/data/fruit-A/"$i"_2.clean.fq.gz;done
# step3 read count
#!bin/bash
for i in $(cat ID.txt); do cp /data/wujing/Myrica/mapping/${i}/${i}Aligned.sortedByCoord.out.bam ./;done
featureCounts -a /data/wujing/Myrica/data/ref/Myrica_rubra_V2.gtf -T 36 -p -g gene_id -t exon -o /data/wujing/Myrica/featureCounts/all_counts.txt /data/wujing/Myrica/bam2/*.sortedByCoord.out.bam

