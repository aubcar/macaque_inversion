#Computational Identification of Structural Variants in Macaque 

###Goals
Use whole genome sequences to computationally identify and analyse structural varients in two species of hybridising macaques

####Samples:   
1. _M. rhesus_ #1:  SRA023856  
    Source: Link to paper in NCBI.   https://www.ncbi.nlm.nih.gov/pubmed/22002653
2. _M. rhesus_ #2: SRR278739 
    Source: Link to paper in NCBI.   
3. _M. fascicularis_: SRA023855  
    Source: Link to paper in NCBI. https://www.ncbi.nlm.nih.gov/pubmed/22002653  

>A 5-year-old female CR macaque and a 4-year-old female CE macaque of Vietnamese origin were used in this study. The CR macaque individual was descended from an individual captured from the wild in Yunnan Province. The origins of these two individuals were confirmed by mitochondrial DNA sequencing. Genomic DNA was collected from the peripheral blood cells of these two individuals.

###Programs used
1.Breakdacer (source)



2. Lumpy (source)



3. Delly (source)


4. Hydra (source) 

###Step1-  Extract SRR files using SRA Toolkit. 
>Locate your study,  download the SRA list( if multiple sequences) and write a forloop to download the fastq files to your workstation and unzip them 


```Shell
for x in {1..125}   #125 comes from the wc -l of Table.txt(taken from the run data)
do
    input=$(awk "NR==$x {print \$0}" ~/rhedata2/Table.txt)
    echo "fastq-dump -I $input --split-files --origfmt" >>$x.sh
    echo "gzip -f ${input}_1.fastq ${input}_2.fastq">>$x.sh 
    chmod a+x $x.sh
done
```

                                                                                                                                           
###Step2 - Alignment and sorting using BWA 0.7.12 and Samtols 1.2

>The sequences must be aligned to the reference genome and sorted into the correct orientation. Putting the puzzle pieces in the correct order to later build the big picture .In our case the reference genome is RheMac3 

```Shell
                                            
for x in {1..125} #125 comes from wc -l of Table.txt                                                                                                                                                    
do
        echo "module load bwa/0.7.12">$x.alignment.sh
        echo "module load samtools/1.2">>$x.alignment.sh
    input=$(awk "NR==$x {print \$0}" ./Table.txt)
    echo "cd /scratch/aubcar/rhe-align/" >>$x.alignment.sh
        #echo "gzip ${input}_1.fastq ${input}_2.fastq" >>$x.alignment.sh                                                                                                                                
        echo "bwa mem -M -v 2 -t 4 rheMac3 ${input}_1.fastq.gz ${input}_2.fastq.gz >$input.bam" >>$x.alignment.sh
        echo "samtools view -Sbu $input.bam | samtools sort - -@ 4 $input.sorted" >>$x.alignment.sh
    chmod a+x $x.alignment.sh
done
```

###Step3 - Adding read group library infomation to the bam files using PICARD tools versin 1.79

>Next read groups were added to each sample. Many algorithms in the GATK pipeline need to know that certain reads were sequenced together in a specific lane, as they attempt to compensate for variability from one sequencing run to the next. Some algorithms need to know the data represents many samples, not just one. Without the read groups and sample information assigned GATK has no way of determining vial information, all reads within a group are assumed to come from the same instrument run, therefore sharing the same error model. Allowing GATK tools to treat all read groups with the same SM value.

```Shell
  for x in {1..125} #125 comes from wc -l of Table.txt
do 
    echo "#! /bin/sh 
source /opt/asn/etc/asn-bash-profiles-special/modules.sh">$x.addRG.sh

        echo "module load picard/1.79">>$x.addRG.sh
    input=$(awk "NR==$x {print \$0}" ./Table.txt)
        library=$(awk "NR==$x+1 {print \$4}" ./SraRunTable.txt)
    echo "cd /scratch/aubcar/rhe-align/" >>$x.addRG.sh
        echo "java -Xms2g -Xmx4g -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/AddOrReplaceReadGroups.jar INPUT=$input.sorted.bam RGID=$input RGLB=$library RGPL=Illumina RGPU=run RGSM=BGI-CR-5 OUTP
UT=$input.wRG.bam">>$x.addRG.sh

    chmod a+x $x.addRG.sh
done
```

###Step4 - Mergingfiles using PICARD tools version 1.79
>A subsequent side effect of all the SRR files is multiple bam files; to make analysis easier, samtools was used to merge the individual lanes back into a single file for easier data manipulation.

```Shell
rhesus_sorted=($(ls SRR*.wRG.bam))

num_files=($(ls SRR*.wRG.bam | wc -l | awk '{print $1}'))

corrected=$((num_files-1))

rm tmp_file
touch tmp_file

for i in $(eval echo {0..$corrected})
        do      
                echo "INPUT=${rhesus_sorted[$i]}" >>tmp_file
        done

inputs=`cat tmp_file`
java -Xms2g -Xmx4g -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/MergeSamFiles.jar $inputs OUTPUT=M_Rhesus.bam VALIDATION_STRINGENCY=SILENT
```

###Step5 - Sorting the newly merged bam with PICARD tools version 1.79
>Bam files were sorted with picard tools as per the GATK guidelines for further data analysis .If samtools is used read groups are sometimes lost causing an error in analysis


```Shell
java -Xms2g -Xmx14g -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/BuildBamIndex.jar I=rheMac3.masked.fa  O=rheMac3.masked.fa.fai
```

###Step6 - Marking the Duplicates using Picard tools version 1.79
>marking the duplicates, this was achieved by using picard tools, the purpose of this step is to mitigate as much as possible the effects of the PCR amplification bias that is introduced during library construction. An added benefit gained downstream is the computational processing is reduced subsequent steps. However there are drawbacks, removing duplicates comes at the cost of setting a rather low hard cap, on the dynamic range of measurement you are going to produce.

```Shell
java -Xms2g -Xmx14g -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/MarkDuplicates.jar INPUT= M_Rhesus.sorted.bam OUTPUT= M_Rhesus.sorted.mkdup.bam METRICS_FILE=Rhesus.txt OPTICAL_DUPLICATE_PIXEL_DIS
TANCE=2500 
```

#Step 7 - Step one, marking the indels and creating an intervals file. Using GATK 
> Step one of indel realignment creates a target file locating the indels for step two to go locally rearange them. 
```Shell
java -Xms2g -Xmx14g -jar /opt/asn/apps/gatk_3.4-46/GenomeAnalysisTK.jar -R rheMac3.masked.fa  -T RealignerTargetCreator -I M_Rhesus.sorted.mkdup.bam  -o M_Rhesus.sorted.mkdup.intervals
```

###Step 8 Indel realignment step 2, using GATK 
>insertion and deletion realignment was run via GATK, indel realigners locally rearrange reads inregions where INDELS might be present in order to more easily identify them.
```Shell
java -Xms2g -Xmx14g -jar /opt/asn/apps/gatk_3.4-46/GenomeAnalysisTK.jar -R rheMac3.masked.fa -T IndelRealigner  -I M_Rhesus2.sorted.mkdup.bam  -o M_Rhesus2.wrg.sorted.mkdup.indel.bam  -targetIntervals M_Rhesus2.sorted.mkdup.intervals
```


###STEP 9 BQSR
```Shelll 


###Change into correct directory
wd=/home/car0019/hbqsr/
cd ${wd}


##load in all required modules
module load vcftools/v0.1.14-14
module load gatk/3.6
module load R/3.2.4
module load samtools

###define the paths to tools being used
gatk_path=/tools/gatk-3.6/
vcf_tools_path=/tools/vcftools-v0.1.14-14/bin/
export PERL5LIB=/tools/vcftools-v0.1.14-14/share/perl5/


###define all variables
reference=rheMac3.masked.fa
bam_file=M_Rhesus.wrg.sorted.mkdup.indel.recal1.bam
filtered_output=var.rhesus1
iteration=second
Qual=500
min_depth=50
max_depth=123
pre_GATK_vcf=${filtered_output}${iteration}.vcf
ouputbam=M_Rhesus.wrg.sorted.mkdup.indel.recal2.bam

###Haplotype caller (step1) get raw vcf

java -Xms2g -Xmx4g -jar ${gatk_path}GenomeAnalysisTK.jar -T HaplotypeCaller -R ${wd}${reference} -I ${bam_file} -o ${filtered_output}${iteration}.vcf

###filter the VCF based on calculated stats

${vcf_tools_path}vcftools --vcf  ${filtered_output}${iteration} --minDP ${min_depth} --maxDP ${max_depth} --max-missing 1 --minQ ${Qual} --recode --recode-INFO-all --out ${wd}${pre_GATK_vcf}
                                                                                                                                                                                                                                                                                                                                         


###make pre recal table using GATK 

java -jar ${gatk_path}GenomeAnalysisTK.jar -T BaseRecalibrator -nct 8 -R ${wd}${reference} -I ${wd}${bam_file} -knownSites ${wd}${pre_GATK_vcf}.recode.vcf -o ${wd}${filtered_output}_recal_data\
.table

                                                                                                                                                     
echo "before gatk" > test.txt

###post table 
java -jar ${gatk_path}GenomeAnalysisTK.jar -T BaseRecalibrator -nct 8 -R ${wd}${reference} -I ${wd}${bam_file} -knownSites ${wd}${pre_GATK_vcf}.recode.vcf -BQSR ${wd}${filtered_output}_recal_d\
ata.table -o ${wd}${filtered_output}_post_recal_data.table

echo "after gatk" >> test.txt


echo "before" > test.txt

###make the plots of the tables

java -jar ${gatk_path}GenomeAnalysisTK.jar -T AnalyzeCovariates -R ${wd}${reference}-l DEBUG -csv ${wd}${filtered_output}.csv -before ${wd}${filtered_output}_recal_data.table -after ${wd}${filtered_output}_post_recal_data.table -plots ${wd}${filtered_output}.recalibration.plots.pdf

echo "after" >> test.txt

echo "before" > test.txt
###add the recalibration to the bam file

java -Xms2g -Xmx8g -jar ${gatk_path}GenomeAnalysisTK.jar -T PrintReads -R ${wd}${reference} -I ${wd}${bam_file} -BQSR ${wd}${filtered_output}_recal_data.table -o ${wd}${outputbam}

echo "after" >>test.txt

###vcf compare to see if you need to run BQSR another time

java -Xms2g -Xmx4g -jar ${gatk_path}GenomeAnalysisTK.jar -T HaplotypeCaller -R ${wd}${reference} -I ${outputbam} -o ${filtered_output}.postrecal.vcf

${vcf_tools_path}vcftools --vcf ${filtered_output}.postrecal.vcf --minDP ${min_depth} --maxDP ${max_depth} --max-missing 1 --minQ ${Qual} --recode --recode-INFO-all --out ${wd}${filtered_output}.postrecal.filtered

file1=${pre_GATK_vcf}.recode.vcf
file2=${filtered_output}.postrecal.filtered.recode.vcf
bgzip ${file1}
bgzip ${file2}
tabix ${file1}
tabix ${file2}
${vcf_tools_path}vcf-compare ${file1} ${file2} > vcf_compare2.txt

###Repeat BQSR till VCF compare=99%
```

###Step 10 extraction Using samtools view version 1.2
>To cut down on computational time and frivilous file sizes down the pipeline, chromosomes 2 and 5 were extracted. 

```Shell
samtools index M_Rhesus.sorted.bam
samtools view -bh M_Rhesus.sorted.bam chr5 >M_Rhesus.sorted.chr5.bam
samtools view -bh M_Rhesus.sorted.bam chr2 >M_Rhesus.sorted.chr2.bam

#if [ -e M_Fascicularis.sorted.bam ] && [ -s M_Fascicularis.sorted.bam ] && [ M_Fascicularis.sorted.bam -nt M_Fascicularis.bam ]; then
#    rm M_Fascicularis.bam
#else
#    echo "error with sorting bam"
#    exit
#fi

```


###Step11
HYDRA/LUMPY/DELLY/BD

####hydra
```shell
####step1 aligment hydra

for x in {1..125} #125 comes from wc -l of Table.txt
do

echo "queue info">$x.hydra.sh
        
        echo "module load bwa/0.7.15">>$x.hydra.sh
        echo "module load samtools/1.3.1">>$x.hydra.sh
        input=$(awk "NR==$x {print \$0}" ./Table.txt)
        echo "bwa aln -t 4 /home/car0019/hydra_rhesus/rheMac3 /home/car0019/hydra_rhesus/${input}_1.fastq > /home/car0019/hydra_rhesus/${input}_1.sai" >> $x.hydra.sh
        echo "bwa aln -t 4 /home/car0019/hydra_rhesus/rheMac3 /home/car0019/hydra_rhesus/${input}_2.fastq > /home/car0019/hydra_rhesus/${input}_2.sai" >> $x.hydra.sh
        echo " bwa sampe /home/car0019/hydra_rhesus/rheMac3 /home/car0019/hydra_rhesus/${input}_1.sai /home/car0019/hydra_rhesus/${input}_2.sai /home/car0019/hydra_rhesus/${input}_1.fastq /home/car0019/hydra_rhesus/${input}_2.fastq  | samtools view -Sbu - > /home/car0019/hydra_rhesus/$input.bam" >> $x.hydra.sh
        chmod a+x $x.hydra.sh
done


####step2  Tier 2 alignment
module load samtools 
module load novocraft/3.4.6
module load bedtools/2
ls | grep ".bam$" > list.txt
files=`cat list.txt` 
for i in ${files} 
do 

samtools view -uF 2 ${i}| bamToFastq -i ${i}  -fq ${i}.ter1.disc.1.fq  -fq2 ${i}.tier1.disc.2.fq
novoalign -d rheMac3 -f ${i}.ter1.disc.1.fq ${i}.tier1.disc.2.fq -i 375 1800 -r Random -o SAM | samtools view -Sb - > ${i}.tier2.queryorder.bam 

done

rm ${i}.ter1.disc.1.fq
rm ${i}.tier1.disc.2.fq


###step3 alignment 
```

####lumpy
```shell
data_dir=/home/car0019/asc_data/data/

ls | grep ".bam$" > list.txt
files=`cat list.txt` 
for i in ${files} 
do 
lumpyexpress -B ${data_dir}${i} -v -P 
done
```

####delly 

```shell
module load delly
wd=/home/car0019/delly/
cd ${wd}
data_dir=/home/car0019/asc_data/data/
ls | grep ".bam$" > list.txt
files=`cat list.txt`
ref=rheMac3.masked.fa 
for i in ${files} 
do 
delly call -g ${data_dir}${ref} -t INV -e illumina ${i} -o ${i}.bcf
done 
```
####convert delly bcf to vcf using bcftools
```
ls | grep ".bcf$" > list.txt
files=`cat list.txt`

for i in ${files} 
do 
bcftools convert ${i} -O v -o {i}.vcf
done
```


####breakdancer
```shell
bam2cfg.pl -g -h -v 2 -q 40 M_Rhesus.wrg.sorted.mkdup.indel.bam > M_Rhesus.wrg.sorted.mkdup.indel.cfg
```

```shell
breakdancer-max -d reads M_Rhesus.wrg.sorted.mkdup.indel.cfg > M_Rhesus.wrg.sorted.mkdup.indel.ctx
```

###Step12-filtering + bedfile formation
> Files were filtered for <1Mb, 90% confidence, and inversion variance and piped into a bed file.

```shell

awk '$0!~/^#/ && $8>=100000 && $8<=90000000 && $9>=90 && $7~/INV/ {print $1,$2,$5}' M_Rhesus.wrg.sorted.mkdup.indel.chr5.ctx >M_Rhesus.wrg.sorted.mkdup.indel.chr5.bed
```


###Step13- analysis 
> Files were loaded into UCSC genome browser and IGV/Tablet for analysis. 


