#!/usr/bin/env sh

# 1000G_omni2.5.hg38
wget -P config/reference https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz
wget -P config/reference https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_omni2.5.hg38.vcf.gz.tbi

# 1000G_phase1.snps.high_confidence.hg38
wget -P config/reference https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz
wget -P config/reference https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/1000G_phase1.snps.high_confidence.hg38.vcf.gz.tbi

# Axiom_Exome_Plus.genotypes.all_populations.poly.hg38
wget -P config/reference https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz
wget -P config/reference https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz.tbi

# GATK executable
wget -P config/reference https://github.com/broadinstitute/gatk/releases/download/4.3.0.0/gatk-4.3.0.0.zip
unzip gatk-4.3.0.0.zip && rm gatk-4.3.0.0.zip

# GnomAD_Genomes
wget -P config/reference https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz
wget -P config/reference https://krishna.gs.washington.edu/download/CADD/v1.6/GRCh38/gnomad.genomes.r3.0.indel.tsv.gz.tbi

# GRCh38_Verily_v1.genome.fa
wget -P config/reference https://storage.googleapis.com/genomics-public-data/references/GRCh38_Verily/GRCh38_Verily_v1.genome.fa
wget -P config/reference https://storage.googleapis.com/genomics-public-data/references/GRCh38_Verily/GRCh38_Verily_v1.genome.fa.fai
wget -P config/reference https://storage.googleapis.com/genomics-public-data/references/GRCh38_Verily/bwa-0.7.12/GRCh38_Verily_v1.genome.fa.alt
wget -P config/reference https://storage.googleapis.com/genomics-public-data/references/GRCh38_Verily/bwa-0.7.12/GRCh38_Verily_v1.genome.fa.amb
wget -P config/reference https://storage.googleapis.com/genomics-public-data/references/GRCh38_Verily/bwa-0.7.12/GRCh38_Verily_v1.genome.fa.ann
wget -P config/reference https://storage.googleapis.com/genomics-public-data/references/GRCh38_Verily/bwa-0.7.12/GRCh38_Verily_v1.genome.fa.bwt
wget -P config/reference https://storage.googleapis.com/genomics-public-data/references/GRCh38_Verily/bwa-0.7.12/GRCh38_Verily_v1.genome.fa.pac
wget -P config/reference https://storage.googleapis.com/genomics-public-data/references/GRCh38_Verily/bwa-0.7.12/GRCh38_Verily_v1.genome.fa.sa

# Hapmap_3.3.hg38
wget -P config/reference https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz
wget -P config/reference https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi

# HG002_bed
wget -P config/reference https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed

# Homo_sapiens_assembly38.dbsnp138
wget -P config/reference https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz
wget -P config/reference https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.gz.tbi

# Homo_sapiens_assembly38.known_indels
wget -P config/reference https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz
wget -P config/reference https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi

# Mills_and_1000G_gold_standard.indels.hg38
wget -P config/reference https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
wget -P config/reference https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi

# Indexed VEP cache v107
wget -P config/reference/v107 https://ftp.ensembl.org/pub/release-107/variation/indexed_vep_cache/homo_sapiens_merged_vep_107_GRCh38.tar.gz
wget -P config/reference/v107 https://ftp.ensembl.org/pub/release-107/variation/indexed_vep_cache/homo_sapiens_refseq_vep_107_GRCh38.tar.gz

# Whole_genome_SNVs
wget -P config/reference https://krishna.gs.washington.edu/download/CADD/v1.4/GRCh38/whole_genome_SNVs.tsv.gz
wget -P config/reference https://krishna.gs.washington.edu/download/CADD/v1.4/GRCh38/whole_genome_SNVs.tsv.gz.tbi
