import os
import sys
import tempfile
from pathlib import Path
from .ete import parse_kraken, update_ncbi_taxonomy, write_kraken_summary


def main(reports: list[Path], summary: Path) -> None:
    report_counts = {
        report.stem.replace("-taxa", ""): parse_kraken(report) for report in reports
    }
    with tempfile.TemporaryDirectory() as tmpdir:
        os.environ["ETE_CACHE_DIR"] = tmpdir
        ncbi = update_ncbi_taxonomy()
        write_kraken_summary(report_counts, summary, ncbi)


if "snakemake" in globals():
    log = snakemake.log[0]  # type: ignore
    reports = [Path(fp) for fp in snakemake.input.reports]  # type: ignore
    summary = Path(snakemake.output.summary)  # type: ignore

    with open(log, "w") as log_f:
        try:
            sys.stdout = log_f
            sys.stderr = log_f
            main(reports, summary)
        except Exception as e:
            log_f.write(f"Error: {e}\n")
            raise


if __name__ == "__main__":
    if len(sys.argv) < 3:
        print(
            "Usage: python summarize_kraken2_reports.py <output_summary.tsv> <report1.tsv> [<report2.tsv> ...]"
        )
        sys.exit(1)

    output_summary = Path(sys.argv[1])
    report_files = [Path(fp) for fp in sys.argv[2:]]

    main(report_files, output_summary)
