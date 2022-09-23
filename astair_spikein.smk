"""
Workflow to compute taps conversion(notrim)
inputs: 
    fastq
outputs:
    meth sta
run:
module load snakemake/5.26.1-foss-2019b-Python-3.7.4 
snakemake -np --use-envmodules --max-status-checks-per-second 0.01 -j 16 --snakefile code/astair_spikein_mESC.smk --cluster "qsub -o astair_spikein_mESC.log -e astair_spikein_mESC.err  -P ludwig.prjc -q long.qc -cwd -V -S /bin/bash -N caps -pe shmem 4"
"""


from snakemake.utils import min_version
min_version("5.26")

configfile: "code/configure.yaml"
SAMPLES = config["CAPS_GBM"]
#SAMPLES = config["CAPS_mESC"]


print(SAMPLES)

rule all:
    input:
        expand("fastq/{sample}_fastp_{readDirection}.fq.gz", sample = SAMPLES, readDirection=['1','2']),
        expand("align/{sample}.fastp.bwa.bam", sample = SAMPLES),
        expand("align/{sample}.md.bam", sample = SAMPLES),
        expand("align/{sample}.spikein.sort.bam", sample = SAMPLES),
        expand("astair_output/{sample}.spikein.sort_mCtoT_all.mods.gz", sample = SAMPLES),
        expand("astair_output/{sample}.all.sta", sample = SAMPLES)




rule fastp: 
    input:
        expand("fastq/{{sample}}_R{readDirection}.fastq.gz",readDirection=['1','2'])
    output:
        expand("fastq/{{sample}}_fastp_{readDirection}.fq.gz", readDirection=['1','2'])
    log:
        "logs/{sample}.trim.log"
    params:
        "fastq/{sample}.fastp"
    threads: 3
    shell:
        """
         /users/ludwig/ebu571/.conda/envs/env_fastp/bin/fastp -i {input[0]} \
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
    threads: 6
    shell:
        """
        (bwa mem -t {threads} {params.ref} {input} -I 500,120,1000,20 |\
        samtools sort -@ 8 -O BAM -T {params.tmp} >{output}) 1>{log} 2>&1
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

################################ Extract trimmed reads on spikein #################################
rule bam_spikein:
    input:
        "align/{sample}.fastp.bwa.bam"
    output:
        "align/{sample}.spikein.sort.bam"
    params:
        ref=REF

    envmodules: "samtools/1.8-gcc5.4.0"
    shell:
        """
        samtools view -bS {input} -L <(awk 'BEGIN{{OFS="\\t";FS="\\t"}}{{if($1!~/chr/)print $1,"1",$2}}' {params.ref}.fai) >{output}
        """

rule astair_spikein:
    input:
        "align/{sample}.spikein.sort.bam"  
    output: 
        "astair_output/{sample}.spikein.sort_mCtoT_all.mods.gz"
    log:
        "logs/astaircall/{sample}.spikein.sort_mCtoT_all.mods.log"

    threads: 8
    shell:
        """
        module load Python/3.7.4-GCCcore-8.3.0
        (/gpfs2/well/ludwig/users/ebu571/Python/asTair_p3.7.4-skylake/bin/astair call -sc True \
        -i {input} -f "/users/ludwig/ebu571/ebu571/14Jul2022_2022_100cycles/resource/spikein.fa" \
        -mq 10 \
        -t 8 -m mCtoT --max_depth 1000000\
        --gz -d astair_output/) 2> {log}
        """

################################ statistics #################################
rule meth_sta_fastp_ncnn:
    input:
         "astair_output/{sample}.spikein.sort_mCtoT_all.mods.gz"

    output: 
        cg="astair_output/{sample}.spikein.sort_mCtoT_CpG.mods.gz",
        chh="astair_output/{sample}.spikein.sort_mCtoT_CHH.mods.gz",
        chg="astair_output/{sample}.spikein.sort_mCtoT_CHG.mods.gz",
        sta="astair_output/{sample}.all.sta"
    shell: 
        """
            zcat {input} | awk '$10=="CpG"' | gzip > {output.cg}
            zcat {input} | awk '$10=="CHH"' | gzip > {output.chh}
            zcat {input} | awk '$10=="CHG"' | gzip > {output.chg}

            ntrim_mCG_144hmC=`zcat {output.cg} | awk -F "\\t" '$1=="144hmC" && $5+$6>0' |awk '$2==86||$2==91||$2==99||$2==109'|awk -F "\\t" '{{mC+=$5;uC+=$6}}END{{if(mC+uC>0)printf("%.4f", mC/(mC+uC));else print "NA"}}'`
            ntrim_nCG_144hmC=`zcat {output.cg} | awk -F "\\t" '$1=="144hmC" && $5+$6>0' |awk '$2==86||$2==91||$2==99||$2==109'|awk -F "\\t" '{{mC+=$5;uC+=$6}}END{{if(mC+uC>0)printf("%.4f", (mC+uC)/NR);else print "NA"}}'`
            ntrim_mCAA_144hmC=`zcat {output.chh}| awk -F "\\t" '$1=="144hmC" && $5+$6>0' |awk '$2==116'|awk -F "\\t" '{{mC+=$5;uC+=$6}}END{{if(mC+uC>0)printf("%.4f", mC/(mC+uC));else print "NA"}}'`
            ntrim_nCAA_144hmC=`zcat {output.chh}| awk -F "\\t" '$1=="144hmC" && $5+$6>0' |awk '$2==116'|awk -F "\\t" '{{mC+=$5;uC+=$6}}END{{if(mC+uC>0)printf("%.4f", (mC+uC)/NR);else print "NA"}}'`
            
            ntrim_mCG_lambda=`zcat {output.cg}|awk -F "\\t" '$1=="J02459.1" && $5+$6>0' |awk -F "\\t" '{{mC+=$5;uC+=$6}}END{{printf("%.4f", mC/(mC+uC))}}'`
            ntrim_nCG_lambda=`zcat {output.cg}|awk -F "\\t" '$1=="J02459.1" && $5+$6>0' |awk -F "\\t" '{{mC+=$5;uC+=$6}}END{{printf("%.4f", (mC+uC)/NR)}}'`
            ntrim_mCH_lambda=`zcat {output.chh} {output.chg} |awk -F "\\t" '$1=="J02459.1" && $5+$6>0'|awk -F "\\t" '{{mC+=$5;uC+=$6}}END{{printf("%.4f", mC/(mC+uC))}}'`
            ntrim_nCH_lambda=`zcat {output.chh} {output.chg} |awk -F "\\t" '$1=="J02459.1" && $5+$6>0'|awk -F "\\t" '{{mC+=$5;uC+=$6}}END{{printf("%.4f", (mC+uC)/NR)}}'`

            ntrim_mnonCG_2kb=`zcat {output.cg}| awk -F "\\t" '$1=="unmodified_2kb" && $5+$6>0' |awk -F "\\t" '{{mC+=$5;uC+=$6}}END{{if(mC+uC>0)printf("%.4f", mC/(mC+uC));else print "NA"}}'`
            ntrim_nnonCG_2kb=`zcat {output.cg}| awk -F "\\t" '$1=="unmodified_2kb" && $5+$6>0' |awk -F "\\t" '{{mC+=$5;uC+=$6}}END{{if(mC+uC>0)printf("%.4f", (mC+uC)/NR);else print "NA"}}'`
            
            ntrim_mnonCH_2kb=`zcat {output.chh} {output.chg}| awk -F "\\t" '$1=="unmodified_2kb" && $5+$6>0' |awk -F "\\t" '{{mC+=$5;uC+=$6}}END{{if(mC+uC>0)printf("%.4f", mC/(mC+uC));else print "NA"}}'`
            ntrim_nnonCH_2kb=`zcat {output.chh} {output.chg}| awk -F "\\t" '$1=="unmodified_2kb" && $5+$6>0' |awk -F "\\t" '{{mC+=$5;uC+=$6}}END{{if(mC+uC>0)printf("%.4f", (mC+uC)/NR);else print "NA"}}'`
            
            ntrim_mnonC_2kb=`zcat {output.cg} {output.chh} {output.chg}| awk -F "\\t" '$1=="unmodified_2kb" && $5+$6>0' |awk -F "\\t" '{{mC+=$5;uC+=$6}}END{{if(mC+uC>0)printf("%.4f", mC/(mC+uC));else print "NA"}}'`
            ntrim_nnonC_2kb=`zcat {output.cg} {output.chh} {output.chg}| awk -F "\\t" '$1=="unmodified_2kb" && $5+$6>0' |awk -F "\\t" '{{mC+=$5;uC+=$6}}END{{if(mC+uC>0)printf("%.4f", (mC+uC)/NR);else print "NA"}}'`
            
            ntrim_mCG_250bp_sym=`zcat {output.cg}|awk -F "\\t" '$1=="250bp" && $5+$6>0' |awk '$2==20||$2==33||$2==35||$2==42||$2==21||$2==34||$2==36||$2==43'|awk -F "\\t" '{{mC+=$5;uC+=$6}}END{{if(mC+uC>0)printf("%.4f", mC/(mC+uC));else print "NA"}}'`
            ntrim_nCG_250bp_sym=`zcat {output.cg}|awk -F "\\t" '$1=="250bp" && $5+$6>0' |awk '$2==20||$2==33||$2==35||$2==42||$2==21||$2==34||$2==36||$2==43'|awk -F "\\t" '{{mC+=$5;uC+=$6}}END{{if(mC+uC>0)printf("%.4f", (mC+uC)/NR);else print "NA"}}'`
            ntrim_mCG_250bp_asym=`zcat {output.cg}|awk -F "\\t" '$1=="250bp" && $5+$6>0' |awk '$2==85||$2==87||$2==96||$2==98||$2==108||$2==115||$2==117||$2==119||$2==128||$2==148||$2==158||$2==169||$2==178||$2==193||$2==201||$2==218||$2==226||$2==228||$2==64||$2==28'|awk -F "\\t" '{{mC+=$5;uC+=$6}}END{{if(mC+uC>0)printf("%.4f", mC/(mC+uC));else print "NA"}}'`
            ntrim_nCG_250bp_asym=`zcat {output.cg}|awk -F "\\t" '$1=="250bp" && $5+$6>0' |awk '$2==85||$2==87||$2==96||$2==98||$2==108||$2==115||$2==117||$2==119||$2==128||$2==148||$2==158||$2==169||$2==178||$2==193||$2==201||$2==218||$2==226||$2==228||$2==64||$2==28'|awk -F "\\t" '{{mC+=$5;uC+=$6}}END{{if(mC+uC>0)printf("%.4f", (mC+uC)/NR);else print "NA"}}'`

            echo -e "smp\\tntrim_mCG_144hmC\\tntrim_nCG_144hmC\\tntrim_mCAA_144hmC\\tntrim_nCAA_144hmC\\tntrim_mCG_250bp_sym\\tntrim_nCG_250bp_sym\\tntrim_mCG_250bp_asym\\tntrim_nCG_250bp_asym\\tntrim_mnonCG_2kb\\tntrim_nnonCG_2kb\\tntrim_mnonCH_2kb\\tntrim_nnonCH_2kb\\tntrim_mnonC_2kb\\tntrim_nnonC_2kb\\tntrim_mCG_lambda\\tntrim_nCG_lambda\\tntrim_mCH_lambda\\tntrim_nCH_lambda" >{output.sta}
            echo -e "{output.sta}\\t${{ntrim_mCG_144hmC}}\\t${{ntrim_nCG_144hmC}}\\t${{ntrim_mCAA_144hmC}}\\t${{ntrim_nCAA_144hmC}}\\t${{ntrim_mCG_250bp_sym}}\\t${{ntrim_nCG_250bp_sym}}\\t${{ntrim_mCG_250bp_asym}}\\t${{ntrim_nCG_250bp_asym}}\\t${{ntrim_mnonCG_2kb}}\\t${{ntrim_nnonCG_2kb}}\\t${{ntrim_mnonCH_2kb}}\\t${{ntrim_nnonCH_2kb}}\\t${{ntrim_mnonC_2kb}}\\t${{ntrim_nnonC_2kb}}\\t${{ntrim_mCG_lambda}}\\t${{ntrim_nCG_lambda}}\\t${{ntrim_mCH_lambda}}\\t${{ntrim_nCH_lambda}}"  >>{output.sta}
        """
