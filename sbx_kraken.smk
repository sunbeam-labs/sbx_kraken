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
        "envs/sbx_kraken_env.yml"
    threads: 8
    shell:
        """
        if LC_ALL=C gzip -l {input[0]} | awk 'NR==2 {{exit($2!=0)}}'; then
            echo "100.00\t0\t0\tR\t1\troot" > {output.report} && \
            echo "C\tA\t1\t136|136\t1:102 |:| 1:102" > {output.raw}
        else
            kraken2 --gzip-compressed \
                    --db {params.db} \
                    --report {output.report} \
                    --output {output.raw} \
                    {params.paired_end} {input} \
                    2>&1 | tee {log}
        fi
        """


rule summarize_kraken2_reports:
    input:
        reports=expand(
            CLASSIFY_FP / "kraken" / "{sample}-taxa.tsv", sample=Samples.keys()
        ),
    output:
        summary=CLASSIFY_FP / "kraken" / "all_samples.tsv",
    benchmark:
        BENCHMARK_FP / "summarize_kraken2_reports.tsv"
    log:
        LOG_FP / "summarize_kraken2_reports.log",
    script:
        "scripts/summarize_kraken2_reports.py"
