"""
Workflow to for caps analysis
inputs: 
    fastq
outputs:
    meth sta
run:
module load snakemake/5.26.1-foss-2019b-Python-3.7.4 
snakemake -np --use-envmodules --max-status-checks-per-second 0.01 -j 16 --snakefile code/standard_caps_mESC.smk --cluster "qsub -o caps_merge_mESC.log -e caps_merge_mESC.err  -P ludwig.prjc -q short.qc@@short.hga -cwd -V -S /bin/bash -N caps_merge -pe shmem 4"
"""


from snakemake.utils import min_version
min_version("5.26")

configfile: "code/configure.yaml"

SAMPLES = config["CAPS_mESC"]
REF = config["MM9_spikein_ref"]
Inputbam = "CAPS_mESC_merged_bam"
print(SAMPLES)

rule all:
    input:
        expand("align/{sample}.md.bam", sample = SAMPLES),
        expand("align/{sample}.md.metrics.txt", sample = SAMPLES),
        expand("align/{sample}.md.q10.bam", sample = SAMPLES),
        expand("align/{sample}.md.q10.bam.bai", sample = SAMPLES),
        expand("align/{sample}.md.bam.bai", sample = SAMPLES),
        expand("align/{sample}.mapping_rate_report.txt", sample = SAMPLES),
        expand("astair_output/{sample}.md.q10_mCtoT_all.mods.gz",  sample = SAMPLES),
        expand("astair_output/{sample}.md.q10_mCtoT_all.chr_rmblack.bedGraph",  sample = SAMPLES),
        expand("align/{inputbam}.bam", inputbam=Inputbam),
        expand("align/{inputbam}.md.bam", inputbam=Inputbam),
        expand("align/{inputbam}.md.q10.bam", inputbam=Inputbam),
        expand("align/{inputbam}.md.q10.bam.bai", inputbam=Inputbam),
        expand("align/{inputbam}.md.bam.bai", inputbam=Inputbam),
        expand("align/{inputbam}.mapping_rate_report.txt", inputbam=Inputbam),
        expand("astair_output/{inputbam}.md.q10_mCtoT_all.mods.gz", inputbam=Inputbam),
        expand("astair_output/{inputbam}.md.q10_mCtoT_all.chr_rmblack.bedGraph", inputbam=Inputbam),
        
rule markdup: 
    input:
        "align/{sample}.fastp.bwa.bam"
    output:
        bam="align/{sample}.md.bam",
        matrix="align/{sample}.md.metrics.txt"
    log:
        "logs/{sample}.picard.log"
    params:
        tmp=temp("markdup") 
    envmodules: "picard/2.23.0-Java-11"
    shell:
        """
        (java -jar /apps/eb/skylake/software/picard/2.23.0-Java-11/picard.jar MarkDuplicates \
            I={input} \
            O={output.bam} \
            M={output.matrix} \
            TMP_DIR={params.tmp} ) 1>{log} 2>&1
        """
rule indexbam: 
    input:
        "align/{sample}.md.bam"
    output:
        q10bam="align/{sample}.md.q10.bam",
        q10idx="align/{sample}.md.q10.bam.bai", 
        q0idx="align/{sample}.md.bam.bai", 
    envmodules: "samtools/1.8-gcc5.4.0"
    params:
        mq=" -q 10 "
    shell:
        """
        samtools view -bS {params.mq} {input} >{output.q10bam}
        samtools index {output.q10bam}
        samtools index {input}
        """

###################### Mapping statistics ###############
rule qc_mapped_reads:
    input:
        fq1="fastq/{sample}_R1.fastq.gz",
        fq2="fastq/{sample}_R2.fastq.gz",
        tm1="fastq/{sample}_fastp_1.fq.gz",
        tm2="fastq/{sample}_fastp_2.fq.gz",
        mapped="align/{sample}.fastp.bwa.bam",
        dedup="align/{sample}.md.bam"
    output:
        "align/{sample}.mapping_rate_report.txt"        
    
    envmodules: "samtools/1.8-gcc5.4.0"
    
    shell:
        """
        echo 'raw reads' > {output}
        echo $(zcat {input.fq1}|wc -l)/4|bc >> {output}
        echo $(zcat {input.fq2}|wc -l)/4|bc >> {output}
        echo 'all reads in fastq' >> {output}
        echo 'trimmed reads' >> {output}
        echo $(zcat {input.tm1}|wc -l)/4|bc >> {output}
        echo $(zcat {input.tm2}|wc -l)/4|bc >> {output}
        echo 'all reads in trimmed fq' >> {output}

        samtools view {input.mapped} | cut -f 1 | sort | uniq | wc -l >> {output}
        echo 'mapped reads (MAPQ 10)' >> {output}
        samtools view -h -F 0x4 -q 10 {input.mapped} | samtools view -F 0x8 | cut -f 1 | sort | uniq | wc -l >> {output}
        echo 'mapped and deduplicated reads (MAPQ 10)' >> {output}
        samtools view -h -F 0x4 -q 10 {input.dedup} | samtools view -F 0x8 | cut -f 1 | sort | uniq | wc -l >> {output}

        echo 'mapped reads (MAPQ 1)' >> {output}
        samtools view -h -F 0x4 -q 1 {input.mapped} | samtools view -F 0x8 | cut -f 1 | sort | uniq | wc -l >> {output}
        echo 'mapped and deduplicated reads (MAPQ 1)' >> {output}
        samtools view -h -F 0x4 -q 1 {input.dedup} | samtools view -F 0x8 | cut -f 1 | sort | uniq | wc -l >> {output}

        echo 'mapped reads (MAPQ 0)' >> {output}
        samtools view -h -F 0x4 {input.mapped} | samtools view -F 0x8 | cut -f 1 | sort | uniq | wc -l >> {output}
        echo 'mapped and deduplicated reads (MAPQ 0)' >> {output}
        samtools view -h -F 0x4 {input.dedup} | samtools view -F 0x8 | cut -f 1 | sort | uniq | wc -l >> {output}
        
        echo 'PCR and optical duplicates(MAPQ>10)' >> {output}
        samtools view -h -F 0x4 -q 10 {input.dedup} | samtools view -h -F 0x8 | samtools view -f 1024 | cut -f 1 | sort | uniq | wc -l >> {output}

        echo 'properly mapped reads (MAPQ 10)' >> {output}
        samtools view -h -F 0x4 -q 10 {input.mapped} | samtools view -h -F 0x8 | awk '$2==83 || $2==99'| cut -f 1 | sort | uniq | wc -l >> {output}

        """    
######################

rule astair_call:
    input: 
        bam="align/{sample}.md.q10.bam",
        bai="align/{sample}.md.q10.bam.bai"
    output:
        "astair_output/{sample}.md.q10_mCtoT_all.mods.gz"

    log:
        "logs/astaircall/{sample}md.q10_mCtoT_all.mods.log"

    threads: 8
    shell:
        """
        module load Python/3.7.4-GCCcore-8.3.0
        (/gpfs2/well/ludwig/users/ebu571/Python/asTair_p3.7.4-skylake/bin/astair call -sc True \
        -i {input.bam} --known_snp /users/ludwig/ebu571/ebu571/CAPS_pub/resource/SNV.vcf.gz -f /users/ludwig/ebu571/ebu571/CAPS_pub/resource/mm9_nochrM.fa \
        -mq 10 \
        -t 8 -m mCtoT  \
        --gz -d astair_output/) 2> {log}
        """

rule astair_filter:
    input:
        inall="astair_output/{sample}.md.q10_mCtoT_all.mods.gz"
        
    output:
        outchr="astair_output/{sample}.md.q10_mCtoT_all.chr.bedGraph.gz",
        rmbl="astair_output/{sample}.md.q10_mCtoT_all.chr_rmblack.bedGraph",
        outcg="astair_output/{sample}.md.q10_mCtoT_CpG.chr_rmblack.bedGraph.gz",
        outchg="astair_output/{sample}.md.q10_mCtoT_CHG.chr_rmblack.bedGraph.gz",
        outchh="astair_output/{sample}.md.q10_mCtoT_CHH.chr_rmblack.bedGraph.gz"

    envmodules: "BEDTools/2.30.0-GCC-10.2.0"

    threads: 8
    shell:
        """
        zcat {input} | grep ^chr | sort -k1,1 -k2,2n | gzip > {output.outchr}
        bedtools intersect -v -a {output.outchr} -b resource/mm9-blacklist.bed > {output.rmbl}
        cat {output.rmbl} | awk -F "\t" '$10=="CpG"' | sort -k1,1 -k2,2n | gzip > {outcg}
        cat {output.rmbl} | awk -F "\t" '$10=="CHG"' | sort -k1,1 -k2,2n | gzip > {outchg}
        cat {output.rmbl} | awk -F "\t" '$10=="CHH"' | sort -k1,1 -k2,2n | gzip > {outchh}

        """

rule astair_rawsignal:
    input: 
        "astair_output/{sample}.md.q10_mCtoT_CpG.chr_rmblack.bedGraph.gz"
    output:
        "astair_output/{sample}.md.q10_mCtoT_CpG.chr_rmblack_raw_signal.bed"

    envmodules: "BEDTools/2.30.0-GCC-10.2.0"

    shell:
        """
        cat {input} | awk '$5+$6>0 && $11=="No"' | sort -k1,1 -k2,2n |\
        bedtools map -a resource/mm9_10kb.bed -b - -c 5,6 -o sum -null "NA" > {output}
        """

################ merged ########
rule merge_bams:
    input: 
        expand("align/{sample}.fastp.bwa.bam", sample = SAMPLES)
    output: 
        expand("align/{inputbam}.bam", inputbam=Inputbam)
    params: 
        " -I ".join(expand("align/{sample}.fastp.bwa.bam", sample = SAMPLES))
    envmodules: 
        "GATK/4.1.7.0-GCCcore-8.3.0-Java-11"
    threads: 1
    shell: 
        """
            gatk MergeSamFiles -I {params} -O {output}
        """
            
rule merge_markdup: 
    input:
        expand("align/{inputbam}.bam", inputbam=Inputbam)
    output:
        bam=expand("align/{inputbam}.md.bam", inputbam=Inputbam),
        matrix=expand("align/{inputbam}.md.metrics.txt", inputbam=Inputbam)
    log:
        expand("logs/{inputbam}.picard.log", inputbam=Inputbam)
    params:
        tmp=temp("markdup") 
    envmodules: "picard/2.23.0-Java-11"
    shell:
        """
        (java -jar /apps/eb/skylake/software/picard/2.23.0-Java-11/picard.jar MarkDuplicates \
            I={input} \
            O={output.bam} \
            M={output.matrix} \
            TMP_DIR={params.tmp} ) 1>{log} 2>&1
        """
rule merge_indexbam: 
    input:
        expand("align/{inputbam}.md.bam", inputbam=Inputbam)
    output:
        q10bam=expand("align/{inputbam}.md.q10.bam", inputbam=Inputbam),
        q10idx=expand("align/{inputbam}.md.q10.bam.bai", inputbam=Inputbam),
        q0idx=expand("align/{inputbam}.md.bam.bai", inputbam=Inputbam)
    envmodules: "samtools/1.8-gcc5.4.0"
    params:
        mq=" -q 10 "
    shell:
        """
        samtools view -bS {params.mq} {input} >{output.q10bam}
        samtools index {output.q10bam}
        samtools index {input}
        """

rule merge_astair_call:
    input: 
        bam=expand("align/{inputbam}.md.q10.bam", inputbam=Inputbam),
        bai=expand("align/{inputbam}.md.q10.bam.bai", inputbam=Inputbam)
    output:
        expand("astair_output/{inputbam}.md.q10_mCtoT_all.mods.gz", inputbam=Inputbam)

    log:
        expand("logs/astaircall/{inputbam}md.q10_mCtoT_all.mods.log", inputbam=Inputbam)

    threads: 8
    shell:
        """
        module load Python/3.7.4-GCCcore-8.3.0
        (/gpfs2/well/ludwig/users/ebu571/Python/asTair_p3.7.4-skylake/bin/astair call -sc True \
        -i {input.bam} --known_snp /users/ludwig/ebu571/ebu571/CAPS_pub/resource/SNV.vcf.gz -f /users/ludwig/ebu571/ebu571/CAPS_pub/resource/mm9_nochrM.fa \
        -mq 10 \
        -t 8 -m mCtoT  \
        --gz -d astair_output/) 2> {log}
        """

rule merge_astair_filter:
    input:
        inall=expand("astair_output/{inputbam}.md.q10_mCtoT_all.mods.gz", inputbam=Inputbam)
        
    output:
        outchr=expand("astair_output/{inputbam}.md.q10_mCtoT_all.chr.bedGraph.gz", inputbam=Inputbam),
        rmbl=expand("astair_output/{inputbam}.md.q10_mCtoT_all.chr_rmblack.bedGraph", inputbam=Inputbam),
        outcg="astair_output/{inputbam}.md.q10_mCtoT_CpG.chr_rmblack.bedGraph.gz",
        outchg="astair_output/{inputbam}.md.q10_mCtoT_CHG.chr_rmblack.bedGraph.gz",
        outchh="astair_output/{inputbam}.md.q10_mCtoT_CHH.chr_rmblack.bedGraph.gz"


    envmodules: "BEDTools/2.30.0-GCC-10.2.0"

    threads: 8
    shell:
        """
        zcat {input} | grep ^chr | sort -k1,1 -k2,2n | gzip > {output.outchr}
        bedtools intersect -v -a {output.outchr} -b resource/mm9-blacklist.bed > {output.rmbl}
        cat {output.rmbl} | awk -F "\t" '$10=="CpG"' | sort -k1,1 -k2,2n | gzip > {outcg}
        cat {output.rmbl} | awk -F "\t" '$10=="CHG"' | sort -k1,1 -k2,2n | gzip > {outchg}
        cat {output.rmbl} | awk -F "\t" '$10=="CHH"' | sort -k1,1 -k2,2n | gzip > {outchh}
        """
