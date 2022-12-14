# -*- mode: Snakemake -*-
#
# Rules for running Kraken

TARGET_CLASSIFY = [CLASSIFY_FP / "kraken" / "all_samples.tsv"]


rule all_classify:
    input:
        TARGET_CLASSIFY,


rule kraken2_classify_report:
    input:
        expand(QC_FP / "decontam" / "{{sample}}_{rp}.fastq.gz", rp=Pairs),
    output:
        raw=CLASSIFY_FP / "kraken" / "raw" / "{sample}-raw.tsv",
        report=CLASSIFY_FP / "kraken" / "{sample}-taxa.tsv",
    benchmark:
        BENCHMARK_FP / "kraken2_classify_report_{sample}.tsv"
    params:
        db=Cfg["sbx_kraken"]["kraken_db_fp"],
        paired_end="--paired" if Cfg["all"]["paired_end"] else "",
    threads: 8
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
        expand(CLASSIFY_FP / "kraken" / "{sample}-taxa.tsv", sample=Samples.keys()),
    output:
        CLASSIFY_FP / "kraken" / "all_samples.biom",
    benchmark:
        BENCHMARK_FP / "kraken2_biom.tsv"
    conda:
        "sbx_kraken.yml"
    shell:
        """
        kraken-biom --max D -o {output} {input}
        """


rule classic_k2_biom:
    input:
        CLASSIFY_FP / "kraken" / "all_samples.biom",
    output:
        CLASSIFY_FP / "kraken" / "all_samples.tsv",
    benchmark:
        BENCHMARK_FP / "classic_k2_biom.tsv"
    conda:
        "sbx_kraken.yml"
    shell:
        """
        biom convert -i {input} -o {output} \
        --to-tsv --header-key=taxonomy --process-obs-metadata=taxonomy \
        --output-metadata-id="Consensus Lineage"
        """
