#!/bin/bash

if [ $# -ne 1 ];then
	  echo "#usage: sh $0 [samplename]"
	    exit
    fi

    SAMPLE=$1
    BWA="/data/etc/bwa/bwa"
    SAMTOOLS="/data/etc/samtools/bin/samtools"
    REFERENCE="/data/reference/ucsc.hg19.fasta"
    JAVA="/usr/bin/java"
    PICARD="/data/etc/picard/picard.jar"
    GATK="/data/etc/gatk/GenomeAnalysisTK.jar"
    SNPEFF="/data/etc/snpEff/snpEff.jar"
    MILLS="/data/etc/bundle/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
    A1KG="/data/etc/bundle/1000G_phase1.indels.hg19.sites.vcf"
    DBSNP138="/data/etc/bundle/dbsnp_138.hg19.vcf"
