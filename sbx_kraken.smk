try:
    SBX_KRAKEN_VERSION = get_ext_version("sbx_kraken")
except NameError:
    # For backwards compatibility with older versions of Sunbeam
    SBX_KRAKEN_VERSION = "0.0.0"
CLASSIFY_FP = output_subdir(Cfg, "classify")


localrules:
    all_kraken,


rule all_kraken:
    input:
        CLASSIFY_FP / "kraken" / "all_samples.tsv",
        CLASSIFY_FP / "bracken" / "all_samples_bracken.tsv",

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
    container:
        f"docker://sunbeamlabs/sbx_kraken:{SBX_KRAKEN_VERSION}"
    threads: 8
    shell:
        """
        if LC_ALL=C gzip -l {input[0]} | awk 'NR==2 {{exit($2!=0)}}'; then
            echo "Empty input, filling with dummy data" > {log}
            echo "100.00\t0\t0\tR\t1\troot" > {output.report}
            echo "C\tA\t1\t136|136\t1:102 |:| 1:102" > {output.raw}
        else
            kraken2 \
            --gzip-compressed \
            --db {params.db} \
            --report {output.report} \
            --output {output.raw} \
            {params.paired_end} {input} \
            > {log} 2>&1
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


rule summarize_bracken_reports:
    input:
        reports=expand(
	    CLASSIFY_FP / "bracken" / "{sample}-species-k2raw.tsv", sample=Samples.keys()
        ),
    output:
        summary=CLASSIFY_FP / "bracken" / "all_samples_bracken.tsv",
    benchmark:
        BENCHMARK_FP / "summarize_bracken_reports.tsv"
    log:
        LOG_FP / "summarize_bracken_reports.log",
    script:
        "scripts/summarize_kraken2_reports.py"


rule bracken_species:
    input:
        CLASSIFY_FP / "kraken" / "{sample}-taxa.tsv",
    output: 
        breport=CLASSIFY_FP / "bracken" / "{sample}-species.tsv",    
        kreport=CLASSIFY_FP / "bracken" / "{sample}-species-k2raw.tsv",    
    log:
        LOG_FP / "bracken_species_{sample}.log",
    params:
        db=Cfg["sbx_kraken"]["kraken_db_fp"],
        readlen=Cfg["sbx_kraken"]["read_length_bp"],
    conda:
        "envs/sbx_kraken_env.yml"
    container:
        f"docker://sunbeamlabs/sbx_kraken:{SBX_KRAKEN_VERSION}"
    threads: 1
    shell:
        """
        bracken \
        -i {input} \
        -d {params.db} \
        -o {output.breport} \
        -w {output.kreport} \
        -r {params.readlen} \
        -l S > {log} 2>&1
        """
    

