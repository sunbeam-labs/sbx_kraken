from pathlib import Path
from scripts.summarize_kraken2_reports_f import (
    parse_kraken2_tsv_report,
    create_kraken2_tsv_report,
    write_kraken2_tsv_summary,
)

report_names = [Path(x).stem for x in snakemake.input.reports]
report_counts, consensus_lineages = create_kraken2_tsv_report(
    [parse_kraken2_tsv_report(open(x)) for x in snakemake.input.reports], report_names
)
write_kraken2_tsv_summary(
    report_counts, consensus_lineages, open(snakemake.output.summary, "w")
)
