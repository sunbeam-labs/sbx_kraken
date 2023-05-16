# -*- mode: Snakemake -*-
#
# Rules for running Kraken

TARGET_CLASSIFY = [CLASSIFY_FP / "kraken" / "all_samples.tsv"]


try:
    BENCHMARK_FP
except NameError:
    BENCHMARK_FP = output_subdir(Cfg, "benchmarks")
try:
    LOG_FP
except NameError:
    LOG_FP = output_subdir(Cfg, "logs")


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
    log:
        LOG_FP / "kraken2_classify_report_{sample}.log",
    params:
        db=Cfg["sbx_kraken"]["kraken_db_fp"],
        paired_end="--paired" if Cfg["all"]["paired_end"] else "",
    conda:
        "sbx_kraken_env.yml"
    threads: 8
    shell:
        """
        if LC_ALL=C gzip -l {input} | awk 'NR==2 {exit($2!=0)}'; then
            echo "Empty reads files" > {log}
            echo "0\t0.0\tk__Bacteria; p__; c__; o__; f__; g__; s__" > {output.report}
        else
            kraken2 --gzip-compressed \
                    --db {params.db} \
                    --report {output.report} \
                    --output {output.raw} \
                    {params.paired_end} {input} \
                    2>&1 | tee {log}
        fi
        """


rule kraken2_biom:
    input:
        expand(CLASSIFY_FP / "kraken" / "{sample}-taxa.tsv", sample=Samples.keys()),
    output:
        CLASSIFY_FP / "kraken" / "all_samples.biom",
    benchmark:
        BENCHMARK_FP / "kraken2_biom.tsv"
    log:
        LOG_FP / "kraken2_biom.log",
    conda:
        "sbx_kraken_env.yml"
    shell:
        """
        kraken-biom --max D -o {output} {input} 2>&1 | tee {log}
        """


rule classic_k2_biom:
    input:
        CLASSIFY_FP / "kraken" / "all_samples.biom",
    output:
        CLASSIFY_FP / "kraken" / "all_samples.tsv",
    benchmark:
        BENCHMARK_FP / "classic_k2_biom.tsv"
    log:
        LOG_FP / "classic_k2_biom.log",
    conda:
        "sbx_kraken_env.yml"
    # script:
    #    "biom_to_tsv.py"
    shell:
        """
        biom convert -i {input} -o {output} \
        --to-tsv --header-key=taxonomy --process-obs-metadata=taxonomy \
        --output-metadata-id="Consensus Lineage" \
        2>&1 | tee {log}
        """
