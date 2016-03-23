#Computational Identification of Structural Variants in Macaque 

###Ideology
Use whole genome sequences to computationally identify and analyse structural varients in two species of hybridising macaques

####Samples:   
1. _M. rhesus_ #1:  SRA023856  
    Source: Link to paper in NCBI.   https://www.ncbi.nlm.nih.gov/pubmed/22002653
2. _M. rhesus_ #2: SRR278739 
    Source: Link to paper in NCBI.   
3. _M. fascicularis_: SRA023855  
    Source: Link to paper in NCBI. https://www.ncbi.nlm.nih.gov/pubmed/22002653  

>A 5-year-old female CR macaque and a 4-year-old female CE macaque of Vietnamese origin were used in this study. The CR macaque individual was descended from an individual captured from the wild in Yunnan Province. The origins of these two individuals were confirmed by mitochondrial DNA sequencing. Genomic DNA was collected from the peripheral blood cells of these two individuals.

###Step1. Extract SRR files using SRA Toolkit. 
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
#module load picard/1.79
#java -Xms2g -Xmx4g -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/AddOrReplaceReadGroups.jar INPUT=SRR066239.sorted.bam RGID=SRR066239 RGLB=MACzhDACDFBAPE-44 RGPL=Illumina RGPU=run RGSM=BGI-CR-5 OU
TPUT=SRR066239.wRG.bam  


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
#! /bin/sh 
source /opt/asn/etc/asn-bash-profiles-special/modules.sh

module load picard/1.79

# mv *.fastq.gz /scratch/aubcar/rhe-align

cd /scratch/aubcar/rhe-align 

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
#!/bin/sh
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load picard/1.79
java -Xms2g -Xmx14g -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/BuildBamIndex.jar I=rheMac3.masked.fa  O=rheMac3.masked.fa.fai
```

###Step6 - Marking the Duplicates using Picard tools version 1.79
>marking the duplicates, this was achieved by using picard tools, the purpose of this step is to mitigate as much as possible the effects of the PCR amplification bias that is introduced during library construction. An added benefit gained downstream is the computational processing is reduced subsequent steps. However there are drawbacks, removing duplicates comes at the cost of setting a rather low hard cap, on the dynamic range of measurement you are going to produce.

```Shell
#!/bin/sh
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load picard/1.79

java -Xms2g -Xmx14g -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/MarkDuplicates.jar INPUT= M_Rhesus.sorted.bam OUTPUT= M_Rhesus.sorted.mkdup.bam METRICS_FILE=Rhesus.txt OPTICAL_DUPLICATE_PIXEL_DIS
TANCE=2500 
```

#Step 7 - Step one, marking the indels and creating an intervals file. Using GATK 
> Step one of indel realignment creates a target file locating the indels for step two to go locally rearange them. 
```Shell
#! /bin/bash
# script to run GATK
# setup to use gatk
source /opt/asn/etc/asn-bash-profiles/modules.sh
module load gatk 
java -Xms2g -Xmx14g -jar /opt/asn/apps/gatk_3.4-46/GenomeAnalysisTK.jar -R rheMac3.masked.fa  -T RealignerTargetCreator -I M_Rhesus.sorted.mkdup.bam  -o M_Rhesus.sorted.mkdup.intervals
```

###Step 8 Indel realignment step 2, using GATK 
>insertion and deletion realignment was run via GATK, indel realigners locally rearrange reads inregions where INDELS might be present in order to more easily identify them.
```Shell
#! /bin/bash
# script to run GATK
# setup to use gatk
source /opt/asn/etc/asn-bash-profiles/modules.sh
module load gatk 
java -Xms2g -Xmx14g -jar /opt/asn/apps/gatk_3.4-46/GenomeAnalysisTK.jar -R rheMac3.masked.fa -T IndelRealigner  -I M_Rhesus2.sorted.mkdup.bam  -o M_Rhesus2.wrg.sorted.mkdup.indel.bam  -targetIntervals
 M_Rhesus2.sorted.mkdup.intervals
```

###Step 9 extraction Using samtools view version 1.2
>To cut down on computational time and frivilous file sizes down the pipeline, chromosomes 2 and 5 were extracted. 
```Shell
#!/bin/sh

#queue: medium-serial
#number nodes: 1
#mem: 4gb
#clock-time: 90:00:00 

# Your job requested : cput=100:00:00,mem=32gb,neednodes=1:ppn=8,nodes=1:ppn=8,walltime=18:00:00
# Your job used : cput=12:58:35,mem=32769464kb,vmem=33228816kb,walltime=05:33:44
# Your job's parallel cpu utilization : 29%
# Your job's memory utilization (mem) : 97.66%

module load samtools/1.2

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

###Step10-bam2cfg.pl

##Step11-breakdancer
###Step12-filtering
###Step13-bedfile formation
###Step14- analysis 



