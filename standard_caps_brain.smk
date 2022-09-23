"""
Workflow to compute taps conversion(notrim)
inputs: 
    fastq
outputs:
    meth sta
run:
module load snakemake/5.26.1-foss-2019b-Python-3.7.4 
snakemake -np --use-envmodules --max-status-checks-per-second 0.01 -j 12 --snakefile code/standard_caps_brain.smk --cluster "qsub -o hg38_samples.log -e hg38_samples.err  -P ludwig.prjc -q long.qc -cwd -V -S /bin/bash -N hg38_samples -pe shmem 3"
"""


from snakemake.utils import min_version
min_version("5.26")

configfile: "code/configure.yaml"


SAMPLES = config["HUMAN_SMP"]
REF = config["HUMAN_REF"]
print(SAMPLES)

rule all:
    input:
        expand("fastq/{sample}_fastp_{readDirection}.fq.gz", sample = SAMPLES, readDirection=['1','2']),
        expand("align/{sample}.fastp.bwa.bam", sample = SAMPLES),
        expand("align/{sample}.md.bam", sample = SAMPLES),
        expand("meth/{sample}.md_mCtoT_all.mods.gz", sample = SAMPLES),
        expand("meth/{sample}.md.astair_{context}.filtered.bedGraph",sample=SAMPLES, context=['CpG','CHG','CHH']),
        expand("meth/{sample}.updown10k.5k.pdf", sample=SAMPLES),
        expand("align/{sample}.collect_metrics.alignment_summary_metrics", sample = SAMPLES),


rule fastp: 
    input:
        expand("fastq/{{sample}}_R{readDirection}.fastq.gz", readDirection=['1','2'])
    output:
        expand("fastq/{{sample}}_fastp_{readDirection}.fq.gz", readDirection=['1','2'])
    log:
        "logs/{sample}.trim.log"
    params:
        "fastq/{sample}.fastp"
    threads: 3
    shell:
        """
        /gpfs2/well/ludwig/users/cfo155/miniconda2/bin/fastp -i {input[0]} \
            -I {input[1]} \
            -o {output[0]} \
            -O {output[1]} \
            -j {params}.json \
            -h {params}.html
        """

rule align_fastp: 
    input:
        expand("fastq/{{sample}}_fastp_{readDirection}.fq.gz", readDirection=['1','2'])
    output:
        "align/{sample}.fastp.bwa.bam"
    log:
        "logs/{sample}.fastp.bwa.log"
    params:
        ref=REF,
        tmp="{sample}"
    envmodules: "BWA/0.7.17-GCC-8.3.0", "samtools/1.8-gcc5.4.0"
    threads: 5
    shell:
        """
        (bwa mem -t {threads} {params.ref} {input} -I 500,120,1000,20 |\
        samtools sort -@ 8 -O BAM -T {params.tmp} -m 1G --threads 3 >{output}) 1>{log} 2>&1
        """

rule mark_dup:
    input:
        "align/{sample}.fastp.bwa.bam"
    output:
        mdbam="align/{sample}.md.bam",
        matrix="align/{sample}.md.matrix.txt"
    params:
        prefix="align/{sample}"
    envmodules: 
        "picard/2.23.0-Java-11"
    shell:
        """
        java -jar /apps/eb/skylake/software/picard/2.23.0-Java-11/picard.jar MarkDuplicates \
            I={input} \
            O={output.mdbam} \
            M={output.matrix}
        """

rule collect_metrics:
    input:
        "align/{sample}.fastp.bwa.bam"
    output:
        "align/{sample}.collect_metrics.alignment_summary_metrics"
    params:
        prefix="align/{sample}.collect_metrics",
        ref=REF
    envmodules: 
        "picard/2.23.0-Java-11"
    shell:
        """
        java -jar /apps/eb/skylake/software/picard/2.23.0-Java-11/picard.jar CollectMultipleMetrics \
            I={input}  \
            O={params.prefix} \
            R={params.ref}
        """


rule astair_filter_c:
    input:
        "meth/{sample}.md_mCtoT_all.mods.gz"
    output:
        expand("meth/{{sample}}.md.astair_{context}.filtered.bedGraph", context=['CpG','CHG','CHH'])
    envmodules: 
        "BEDTools/2.27.1-foss-2018b"
    params:
        excludregion="resource/hg38_exclude.bed"
    shell: 
        """
        zcat {input}|grep ^chr|awk '$10=="CpG"&&$11=="No"' |awk 'BEGIN{{OFS="\\t"}}{{print $1,$2,$3,$4,$5,$5+$6}}'|sort -k1,1 -k2,2n|\
        bedtools intersect -a - -b {params.excludregion} -v  -sorted >{output[0]}
        zcat {input}|grep ^chr|awk '$10=="CHG"&&$11=="No"' |awk 'BEGIN{{OFS="\\t"}}{{print $1,$2,$3,$4,$5,$5+$6}}'|sort -k1,1 -k2,2n|\
        bedtools intersect -a - -b {params.excludregion} -v  -sorted >{output[1]}
        zcat {input}|grep ^chr|awk '$10=="CHH"&&$11=="No"' |awk 'BEGIN{{OFS="\\t"}}{{print $1,$2,$3,$4,$5,$5+$6}}'|sort -k1,1 -k2,2n|\
        bedtools intersect -a - -b {params.excludregion} -v  -sorted >{output[2]}
        """

rule astair_deeptools:
    input:
        expand("meth/{{sample}}.md.astair_{context}.filtered.bedGraph", context=['CpG','CHG','CHH'])
    output:
        bg=expand("meth/{{sample}}.md.astair_{context}.filtered.chr.temp.bedGraph", context=['CpG','CHG','CHH']),
        bw=expand("meth/{{sample}}.md.astair_{context}.filtered.chr.bw", context=['CpG','CHG','CHH']),
        mat="meth/{sample}.updown10k.5k.mat.gz",
        pdf="meth/{sample}.updown10k.5k.pdf"
    envmodules: 
        "deepTools/3.3.1-foss-2018b-Python-3.6.6", "plotly.py/4.4.1-foss-2018b-Python-3.6.6"
    params:
        region="resource/MANE.GRCh38.v1.0.refseq.gene.bed",
        gsize="resource/GRCh38.gsize.txt"
    shell: 
        """
        grep ^chr[1-9,X,Y] {input[0]} |awk '$6>=5'|cut -f1-4 >{output.bg[0]}
        grep ^chr[1-9,X,Y] {input[1]} |awk '$6>=5'|cut -f1-4 >{output.bg[1]}
        grep ^chr[1-9,X,Y] {input[2]} |awk '$6>=5'|cut -f1-4 >{output.bg[2]}

        bedGraphToBigWig {output.bg[0]} {params.gsize} {output.bw[0]}
        bedGraphToBigWig {output.bg[1]} {params.gsize} {output.bw[1]}
        bedGraphToBigWig {output.bg[2]} {params.gsize} {output.bw[2]}

        computeMatrix scale-regions -S {output.bw} -R {params.region} \
                              --beforeRegionStartLength 5000 \
                              --regionBodyLength 5000 \
                              --afterRegionStartLength 5000 \
                              -o {output.mat} --binSize 10 -p 8

        plotProfile -m {output.mat} \
              -out {output.pdf} \
              --perGroup --samplesLabel healthy GBM --legendLocation "center-right" --colors "#67A9CF" "#EF8A62"
        
        """


rule bam_classify:
    input:
        "align/{sample}.fastp.bwa.bam"
    output:
        chrbam="align/{sample}.chr.sort.bam",
        lambam="align/{sample}.lambda.sort.bam"
    shell:
        """
        samtools view -h {input}|awk '$0~/^@/||$3~/chr/'|samtools view -bS - >{output.chrbam}
        samtools view -h {input}|awk '$0~/^@/||$3~/J02459.1/'|samtools view -bS - >{output.lambam}
        """

rule astair_call: 
    input:
        "align/{sample}.md.bam"
    output:
        meth="meth/{sample}.md_mCtoT_all.mods.gz"
    params:
        ref=REF,
        par=" --minimum_base_quality 13 -d meth -mq 10 --gz -t 5"
    shell:
        """
        /users/ludwig/cfo155/miniconda2/envs/env4py3/bin/astair  call -i {input} -f {params.ref} {params.par} 
        """

rule extract_chr: 
    input:
        "align/{sample}.md.bam"
    output:
        chrbam="align/{sample}.chr.md.bam",
    params:
        ref="/gpfs2/well/ludwig/users/cfo155/14Jul2022/resource/GRCh38_chr.fa",
    shell:
        """
        samtools view -bS {input} -L <(awk 'BEGIN{{OFS="\\t"}}{{print $1,"0",$2}}' {params.ref}.fai)  >{output.chrbam}
        """

