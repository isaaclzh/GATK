import pandas as pd
import pathlib
import os

if not pathlib.Path("config/misc/samples.tsv").exists():
    os.system("python config/config.py")

configfile: "config/config.yml"
gatk = config["gatk_exe"]
ref = config["ref"]["path"]
sample_file = config["sample_file"]
samples = pd.read_table(sample_file)['Sample'].values
chroms = pd.read_csv(config["ref"]["idx"], header=None, sep = '\t').iloc[:,0].values.tolist()[:24]

rule all:
    input:
        expand("results/{sample}/gvcf/{sample}.g.vcf.gz", sample = samples),
        expand("results/annotation/chr{chrom}_recalibrated.vep", chrom = config["ref"]["chr"])

rule bwa_alignment:
    output:
        bam = protected("results/{sample}/bam/{sample}.bwa.bam"),
        index = protected("results/{sample}/bam/{sample}.bwa.bam.bai")
    input:
        forward_reads = expand("fastq/{{sample}}_1{ext}", ext = config["ext"]),
        reverse_reads = expand("fastq/{{sample}}_2{ext}", ext = config["ext"])
    params:
        # Sequencer: Illumina
        rg = r"@RG\tID:{sample}\tSM:{sample}\tPL:Illumina\tLB:{sample}"
    log:
        "logs/bwa/{sample}.log"
    threads:
        config["params"]["bwa"]["threads"]
    shell:
        """
        (bwa mem -M -t {threads} -R '{params.rg}' \
        {ref} {input.forward_reads} {input.reverse_reads} \
        | samtools view -Sb - \
        | samtools sort - -@ {threads} -m 4G -o {output.bam}) 2> {log}
        samtools index {output.bam}
        """

rule mark_duplicates:
    output:
        bam = temp("results/{sample}/bam/{sample}.bwa.sortdup.bam"),
        index = temp("results/{sample}/bam/{sample}.bwa.sortdup.bam.bai"),
        qc = "results/{sample}/qc/duplicates.txt"
    input:
        bam = "results/{sample}/bam/{sample}.bwa.bam",
        index = "results/{sample}/bam/{sample}.bwa.bam.bai"
    params:
        *config["params"]["picard"]["MarkDuplicates"]
    log:
        "logs/picard/sortdup/{sample}.log"
    shell:
        """
        ({gatk} --java-options '-Xmx20G -XX:ParallelGCThreads=8' MarkDuplicates \
        -I {input.bam} -O {output.bam} -M {output.qc} {params}) 2> {log}
        samtools index {output.bam}
        """

rule quality_check:
    output:
        coverage = "results/{sample}/qc/coverage.txt",
        flagstat = "results/{sample}/qc/flagstat.txt"
    input:
        "results/{sample}/bam/{sample}.bwa.sortdup.bam"
    log:
        "logs/picard/collectwgsmetrics/{sample}.log"
    shell:
        """
        ({gatk} --java-options '-Xmx16G' CollectWgsMetrics \
        -I {input} -O {output.coverage} -R {ref}) 2> {log}
        samtools flagstat {input} > {output.flagstat}
        """

rule collate_seq_depth:
    output:
        "config/misc/depth.txt"
    input:
        expand("results/{sample}/qc/coverage.txt", sample = samples)
    script:
        "config/depth.py"

rule recalibrate_base_qualities:
    output:
        temp("results/{sample}/bam/{sample}.recal")
    input:
        bam = "results/{sample}/bam/{sample}.bwa.sortdup.bam",
        index = "results/{sample}/bam/{sample}.bwa.sortdup.bam.bai",
        coverage = "results/{sample}/qc/coverage.txt",
        flagstat = "results/{sample}/qc/flagstat.txt",
        depth = "config/misc/depth.txt"
    params:
        *config["params"]["gatk"]["BaseRecalibrator"]
    log:
        "logs/gatk/bqsr/{sample}.log"
    shell:
        """
        ({gatk} --java-options '-Xmx4G' BaseRecalibrator \
        -I {input.bam} -O {output} -R {ref} {params}) 2> {log}
        """

rule apply_base_quality_recalibration:
    output:
        bqsr = temp("results/{sample}/bam/{sample}.bwa.sortdup.bqsr.bam"),
        index = temp("results/{sample}/bam/{sample}.bwa.sortdup.bqsr.bam.bai")
    input:
        bam = "results/{sample}/bam/{sample}.bwa.sortdup.bam",
        recal = "results/{sample}/bam/{sample}.recal"
    params:
        *config["params"]["gatk"]["ApplyBQSR"]
    log:
        "logs/gatk/apply-bqsr/{sample}.log"
    shell:
        """
        ({gatk} --java-options '-Xmx4G' ApplyBQSR \
        -I {input.bam} -O {output.bqsr} -R {ref} -bqsr {input.recal} {params}) 2> {log}
        samtools index {output.bqsr}
        """

if list(map(lambda chrom:int(chrom) if chrom.isdigit() else chrom, [chrom[3:] for chrom in chroms])) == config["ref"]["chr"]:
    tmp = pd.read_csv(config["ref"]["idx"], header=None, sep = '\t')
    tmp['start'] = 1
    for n, i in enumerate(chroms[::2]):
        partner = chroms[1::2][11-n]
        bed = tmp.loc[tmp.iloc[:,0].isin([i, partner])].loc[:,[0,'start', 1]]
        bed.to_csv('config/misc/%s.bed' % i, sep = '\t', index=False, header=False)

    rule split_bams:
        output:
            bam = temp("results/{sample}/bam/chr{chrom}.bam"),
            index = temp("results/{sample}/bam/chr{chrom}.bam.bai")
        input:
            bqsr = "results/{sample}/bam/{sample}.bwa.sortdup.bqsr.bam",
            index = "results/{sample}/bam/{sample}.bwa.sortdup.bqsr.bam.bai"
        params:
            chrom = config["ref"]["chr"][::2]
        shell:
            """
            samtools view -bh {input.bqsr} --region-file config/misc/chr{wildcards.chrom}.bed > {output.bam}
            samtools index {output.bam}
            """

    rule call_variants:
        output:
            gvcf = temp("results/{sample}/gvcf/chr{chrom}.g.vcf.gz"),
            index = temp("results/{sample}/gvcf/chr{chrom}.g.vcf.gz.tbi")
        input:
            bam = "results/{sample}/bam/chr{chrom}.bam",
            index = "results/{sample}/bam/chr{chrom}.bam.bai"
        params:
            *config["params"]["gatk"]["HaplotypeCaller"]
        log:
            "logs/gatk/haplotypecaller/{sample}/chr{chrom}.log"
        shell:
            """
            ({gatk} --java-options "-Xmx16g" HaplotypeCaller \
            -I {input.bam} -O {output.gvcf} -L config/misc/chr{wildcards.chrom}.bed -R {ref} {params}) 2> {log}
            """

    rule combine_gvcfs:
        output:
            gvcf = protected("results/{sample}/gvcf/{sample}.g.vcf.gz"),
            index = protected("results/{sample}/gvcf/{sample}.g.vcf.gz.tbi")
        input:
            gvcf = expand("results/{{sample}}/gvcf/chr{chrom}.g.vcf.gz", chrom = config['ref']['chr'][::2]),
            index = expand("results/{{sample}}/gvcf/chr{chrom}.g.vcf.gz.tbi", chrom = config['ref']['chr'][::2])
        params:
            *config["params"]["gatk"]["CombineGVCFs"]
        log:
            "logs/gatk/combinegvcfs/{sample}.log"
        run:
            gvcfs = list(map("--variant {}".format, input.gvcf))
            shell("({gatk} CombineGVCFs {gvcfs} -O {output.gvcf} -R {ref} {params}) 2> {log}")

else:
    rule split_bams:
        output:
            bam = temp("results/{sample}/bam/chr{chrom}.bam"),
            index = temp("results/{sample}/bam/chr{chrom}.bam.bai")
        input:
            bam = "results/{sample}/bam/{sample}.bwa.sortdup.bqsr.bam",
            index = "results/{sample}/bam/{sample}.bwa.sortdup.bqsr.bam.bai"
        params:
            chrom = config["ref"]["chr"]
        shell:
            """
            samtools view -bh {input.bam} chr{wildcards.chrom} > {output.bam}
            samtools index {output.bam}
            """

    rule call_variants:
        output:
            gvcf = temp("results/{sample}/gvcf/chr{chrom}.g.vcf.gz"),
            index = temp("results/{sample}/gvcf/chr{chrom}.g.vcf.gz.tbi")
        input:
            bam = "results/{sample}/bam/chr{chrom}.bam",
            index = "results/{sample}/bam/chr{chrom}.bam.bai"
        params:
            *config["params"]["gatk"]["HaplotypeCaller"]
        log:
            "logs/gatk/haplotypecaller/{sample}/chr{chrom}.log"
        shell:
            """
            ({gatk} --java-options "-Xmx16g" HaplotypeCaller \
            -I {input.bam} -O {output.gvcf} -L chr{wildcards.chrom} -R {ref} {params}) 2> {log}
            """

    rule combine_gvcfs:
        output:
            gvcf = protected("results/{sample}/gvcf/{sample}.g.vcf.gz"),
            index = protected("results/{sample}/gvcf/{sample}.g.vcf.gz.tbi")
        input:
            gvcf = expand("results/{{sample}}/gvcf/chr{chrom}.g.vcf.gz", chrom = config['ref']['chr']),
            index = expand("results/{{sample}}/gvcf/chr{chrom}.g.vcf.gz.tbi", chrom = config['ref']['chr'])
        params:
            *config["params"]["gatk"]["CombineGVCFs"]
        log:
            "logs/gatk/combinegvcfs/{sample}.log"
        run:
            gvcfs = list(map("--variant {}".format, input.gvcf))
            shell("({gatk} CombineGVCFs {gvcfs} -O {output.gvcf} -R {ref} {params}) 2> {log}")

rule genotype_gvcfs:
    output:
        vcf = protected("results/{sample}/vcf/{sample}_to_annotate.vcf.gz"),
        index = protected("results/{sample}/vcf/{sample}_to_annotate.vcf.gz.tbi")
    input:
        gvcf = "results/{sample}/gvcf/{sample}.g.vcf.gz",
        index = "results/{sample}/gvcf/{sample}.g.vcf.gz.tbi"
    params:
        *config["params"]["gatk"]["GenotypeGVCFs"]
    log:
        "logs/gatk/genotypegvcfs/{sample}.log"
    shell:
        """
        ({gatk} --java-options "-Xmx5G" GenotypeGVCFs \
        -V {input.gvcf} -O {output.vcf} -R {ref} {params}) 2> {log}
        """

rule genomicsdb:
    output:
        temp(directory("config/database_chr{chrom}"))
    input:
        cohort = "config/misc/cohort.txt",
        vcf = expand("results/{sample}/vcf/{sample}_to_annotate.vcf.gz", sample = samples),
        index_vcf = expand("results/{sample}/vcf/{sample}_to_annotate.vcf.gz.tbi", sample = samples),
        gvcf = expand("results/{sample}/gvcf/{sample}.g.vcf.gz", sample = samples),
        index_gvcf = expand("results/{sample}/gvcf/{sample}.g.vcf.gz.tbi", sample = samples)
    params:
        *config["params"]["gatk"]["GenomicsDBImport"]
    log:
        "logs/gatk/genomicsdb/chr{chrom}.log"
    shell:
        """
        ({gatk} --java-options "-Xmx4G" GenomicsDBImport \
        --sample-name-map {input.cohort} --genomicsdb-workspace-path {output} -L chr{wildcards.chrom} {params}) 2> {log}
        """

rule genotype_cohort:
    output:
        vcf = temp("results/cohort/chr{chrom}_cohort.vcf.gz"),
        index = temp("results/cohort/chr{chrom}_cohort.vcf.gz.tbi")
    input:
        "config/database_chr{chrom}"
    params:
        *config["params"]["gatk"]["GenotypeGVCFs"]
    log:
        "logs/gatk/genotypegvcfs/cohort/chr{chrom}.log"
    shell:
        """
        ({gatk} --java-options "-Xmx5G" GenotypeGVCFs \
        -V gendb://{input} -O {output.vcf} -L chr{wildcards.chrom} --only-output-calls-starting-in-intervals true -R {ref} {params}) 2> {log}
        """

rule cohort_vcf:
    output:
        vcf_cohort = protected("results/cohort/raw_variants.vcf.gz"),
        index_cohort = protected("results/cohort/raw_variants.vcf.gz.tbi")
    input:
        vcf = expand("results/cohort/chr{chrom}_cohort.vcf.gz", chrom = config['ref']['chr']),
        index = expand("results/cohort/chr{chrom}_cohort.vcf.gz.tbi", chrom = config['ref']['chr'])
    params:
        *config["params"]["picard"]["MergeVcfs"]
    log:
        "logs/picard/mergevcfs/cohort.log"
    run:
        vcfs = list(map("-I {}".format, input.vcf))
        shell("({gatk} MergeVcfs {vcfs} -O {output.vcf_cohort} -R {ref} {params}) 2> {log}")

rule recalibrate_snps:
    output:
        recal = temp("results/cohort/recalibrate_snp.recal"),
        index = temp("results/cohort/recalibrate_snp.recal.idx"),
        tranches = temp("results/cohort/recalibrate_snp.tranches")
    input:
        "results/cohort/raw_variants.vcf.gz"
    params:
        *config["params"]["gatk"]["VariantRecalibrator_SNP"]
    log:
        "logs/gatk/recalibrate/recalibrate_snp.log"
    shell:
        """
        ({gatk} --java-options "-Xmx16G" VariantRecalibrator \
        -V {input} -O {output.recal} -R {ref} {params}) 2> {log}
        """

rule apply_vqsr_snps:
    output:
        vcf = temp("results/cohort/recalibrated_snps_raw_indels.vcf.gz"),
        index = temp("results/cohort/recalibrated_snps_raw_indels.vcf.gz.tbi")
    input:
        vcf = "results/cohort/raw_variants.vcf.gz",
        recal = "results/cohort/recalibrate_snp.recal",
        index = "results/cohort/recalibrate_snp.recal.idx",
        tranches = "results/cohort/recalibrate_snp.tranches"
    params:
        *config["params"]["gatk"]["ApplyVQSR_SNP"]
    log:
        "logs/gatk/apply-vqsr/apply_vqsr_snp.log"
    shell:
        """
        ({gatk} --java-options "-Xmx16G" ApplyVQSR \
        -V {input.vcf} -O {output.vcf} -R {ref} \
        --tranches-file {input.tranches} \
        --recal-file {input.recal} {params}) 2> {log}
        """

rule recalibrate_indels:
    output:
        recal = temp("results/cohort/recalibrate_indel.recal"),
        index = temp("results/cohort/recalibrate_indel.recal.idx"),
        tranches = temp("results/cohort/recalibrate_indel.tranches")
    input:
        vcf = "results/cohort/recalibrated_snps_raw_indels.vcf.gz",
        index = "results/cohort/recalibrated_snps_raw_indels.vcf.gz.tbi"
    params:
        *config["params"]["gatk"]["VariantRecalibrator_INDEL"]
    log:
        "logs/gatk/recalibrate/recalibrate_indel.log"
    shell:
        """
        ({gatk} --java-options "-Xmx16G" VariantRecalibrator \
        -V {input.vcf} -O {output.recal} -R {ref} {params}) 2> {log}
        """

rule apply_vqsr_indels:
    output:
        vcf = protected("results/cohort/recalibrated_snps_indels.vcf.gz"),
        index = protected("results/cohort/recalibrated_snps_indels.vcf.gz.tbi")
    input:
        vcf = "results/cohort/recalibrated_snps_raw_indels.vcf.gz",
        index_vcf = "results/cohort/recalibrated_snps_raw_indels.vcf.gz.tbi",
        recal = "results/cohort/recalibrate_indel.recal",
        index_recal = "results/cohort/recalibrate_indel.recal.idx",
        tranches = "results/cohort/recalibrate_indel.tranches"
    params:
        *config["params"]["gatk"]["ApplyVQSR_INDEL"]
    log:
        "logs/gatk/apply-vqsr/apply_vqsr_indel.log"
    shell:
        """
        ({gatk} --java-options "-Xmx16G" ApplyVQSR \
        -V {input.vcf} -O {output.vcf} -R {ref} \
        --tranches-file {input.tranches} \
        --recal-file {input.recal} {params}) 2> {log}
        """

if config['bed'].lower() == 'no':
    rule split_vcfs:
        output:
            vcf = temp("results/cohort/chr{chrom}_recalibrated.vcf.gz"),
            index = temp("results/cohort/chr{chrom}_recalibrated.vcf.gz.tbi")
        input:
            vcf = "results/cohort/recalibrated_snps_indels.vcf.gz",
            index = "results/cohort/recalibrated_snps_indels.vcf.gz.tbi"
        shell:
            """
            bcftools view {input.vcf} -Oz -r chr{wildcards.chrom} -o {output.vcf}
            tabix {output.vcf}
            """

elif config['bed'].lower() == 'yes':
    rule apply_bed_regions:
        output:
            vcf = protected("results/cohort/recalibrated_snps_indels_bed.vcf.gz"),
            index = protected("results/cohort/recalibrated_snps_indels_bed.vcf.gz.tbi")
        input:
            vcf = "results/cohort/recalibrated_snps_indels.vcf.gz",
            index = "results/cohort/recalibrated_snps_indels.vcf.gz.tbi"
        params:
            config["ref"]["bed"]
        shell:
            """
            bcftools view -R {params} {input.vcf} -o {output.vcf}
            tabix {output.vcf}
            """

    rule split_vcfs:
        output:
            vcf = temp("results/cohort/chr{chrom}_recalibrated.vcf.gz"),
            index = temp("results/cohort/chr{chrom}_recalibrated.vcf.gz.tbi")
        input:
            vcf = "results/cohort/recalibrated_snps_indels_bed.vcf.gz",
            index = "results/cohort/recalibrated_snps_indels_bed.vcf.gz.tbi"
        shell:
            """
            bcftools view {input.vcf} -Oz -r chr{wildcards.chrom} -o {output.vcf}
            tabix {output.vcf}
            """

else:
    sys.exit()

rule annotate:
    output:
        "results/annotation/chr{chrom}_recalibrated.vep"
    input:
        vcf = "results/cohort/chr{chrom}_recalibrated.vcf.gz",
        index = "results/cohort/chr{chrom}_recalibrated.vcf.gz.tbi"
    params:
        *config["params"]["vep"]["parameters"],
        *config["params"]["vep"]["plugins"]
    shell:
        """
        vep -i {input.vcf} -o {output} --fasta {ref} {params}
        """
