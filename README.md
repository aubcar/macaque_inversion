# macaque_inversion
Chase Rushton Master's Thesis on macaque speciation

###Step1
```

for x in {1..22}   
do
    input=$(awk "NR==$x {print \$0}" ~/rhedata2/Table2.txt)
    echo "fastq-dump -I $input --split-files --origfmt" >>$x.sh
    echo "gzip -f ${input}_1.fastq ${input}_2.fastq">>$x.sh 
    chmod a+x $x.sh
done
```

 nodes 4,8gbs ram,48hrs, medium-parallel                                                                                                                                                               
###Step2
```
#cd ~/                                                                                                                                                                                                  
#cp  stevison/rheMac3.*  /scratch/aubcar/rhe-align/                                                                                                                                                     
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

###Step3
```
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

# 5 hours, 1 cpu , 4gb small-serial 
```

###Step4

```
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


# Ran with options medium-serial, 4gb, 1 cpu , 90 hours 
```

###Step5
```
#!/bin/sh
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load picard/1.70
java -Xms2g -Xmx14g -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/BuildBamIndex.jar I=rheMac3.masked.fa  O=rheMac3.masked.fa.fai
```

###Step6
```
#!/bin/sh
source /opt/asn/etc/asn-bash-profiles-special/modules.sh
module load picard/1.70

java -Xms2g -Xmx14g -jar /opt/asn/apps/picard_1.79/picard-tools-1.79/MarkDuplicates.jar INPUT= M_Rhesus.sorted.bam OUTPUT= M_Rhesus.sorted.mkdup.bam METRICS_FILE=Rhesus.txt OPTICAL_DUPLICATE_PIXEL_DIS
TANCE=2500 
```

#Step 7
```
#! /bin/bash
# script to run GATK
# setup to use gatk
source /opt/asn/etc/asn-bash-profiles/modules.sh
module load gatk 
java -Xms2g -Xmx14g -jar /opt/asn/apps/gatk_3.4-46/GenomeAnalysisTK.jar -R rheMac3.masked.fa  -T RealignerTargetCreator -I M_Rhesus.sorted.mkdup.bam  -o M_Rhesus.sorted.mkdup.intervals
```

###Step 8
```
#! /bin/bash
# script to run GATK
# setup to use gatk
source /opt/asn/etc/asn-bash-profiles/modules.sh
module load gatk 
java -Xms2g -Xmx14g -jar /opt/asn/apps/gatk_3.4-46/GenomeAnalysisTK.jar -R rheMac3.masked.fa -T IndelRealigner  -I M_Rhesus2.sorted.mkdup.bam  -o M_Rhesus2.wrg.sorted.mkdup.indel.bam  -targetIntervals
 M_Rhesus2.sorted.mkdup.intervals
```

###Step 9
```
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

###Step10

