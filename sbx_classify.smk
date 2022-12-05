# -*- mode: Snakemake -*-
#
# Rules for running Kraken

TARGET_CLASSIFY = [CLASSIFY_FP/'kraken'/'all_samples.tsv']

rule all_classify:
    input:
        TARGET_CLASSIFY
        
rule kraken2_classify_report:
    input:
        expand(QC_FP/'decontam'/'{{sample}}_{rp}.fastq.gz', rp = Pairs)
    output:
        raw = CLASSIFY_FP/'kraken'/'raw'/'{sample}-raw.tsv',
        report = CLASSIFY_FP/'kraken'/'{sample}-taxa.tsv'
    params:
        db = Cfg['sbx_classify']['kraken_db_fp'],
        paired_end = "--paired" if Cfg['all']['paired_end'] else ""
    threads:
        Cfg['sbx_classify']['threads']
    shell:
        """
        kraken2 --gzip-compressed \
                --db {params.db} \
                --report {output.report} \
                {params.paired_end} {input} \
                > {output.raw}
        """

rule kraken2_biom:
    input:
        expand(CLASSIFY_FP/'kraken'/'{sample}-taxa.tsv',
               sample=Samples.keys())
    conda:
        "sbx_classify.yml"
    output:
        CLASSIFY_FP/'kraken'/'all_samples.biom'
    shell:
        """
        kraken-biom --max D -o {output} {input}
        """

rule classic_k2_biom:
    input:
        CLASSIFY_FP/'kraken'/'all_samples.biom'
    output:
        CLASSIFY_FP/'kraken'/'all_samples.tsv'
    conda:
        "sbx_classify.yml"
    shell:
        """
        biom convert -i {input} -o {output} \
        --to-tsv --header-key=taxonomy --process-obs-metadata=taxonomy \
        --output-metadata-id="Consensus Lineage"
        """

