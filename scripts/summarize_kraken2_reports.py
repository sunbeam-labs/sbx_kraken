import csv
import os
import sys
import tempfile
from ete3 import NCBITaxa
from pathlib import Path


allowed_ranks = {
    "D": "k__",
    "P": "p__",
    "C": "c__",
    "O": "o__",
    "F": "f__",
    "G": "g__",
    "S": "s__",
}
ncbi_rank_to_prefix = {
    "domain": "k__",
    "phylum": "p__",
    "class": "c__",
    "order": "o__",
    "family": "f__",
    "genus": "g__",
    "species": "s__",
}


def parse_kraken(fp: Path) -> dict[str, str]:
    results = {}
    with open(fp, "r") as f:
        reader = csv.reader(f, delimiter="\t")
        for row in reader:
            (
                percentage,
                fragments_covered,
                fragments_assigned,
                rank,
                taxon_id,
                scientific_name,
            ) = row
            if rank in allowed_ranks:
                results[taxon_id] = fragments_assigned

    return results


def get_consensus_lineage(taxon_id: str, ncbi: NCBITaxa) -> str:
    lineage = ncbi.get_lineage(taxon_id)
    names = ncbi.get_taxid_translator(lineage)
    lineage_ranks = ncbi.get_rank(lineage)

    if not lineage:
        raise ValueError(f"Taxon ID {taxon_id} not found in NCBI taxonomy database.")

    lineage_parts = []
    for taxid in lineage:
        rank = lineage_ranks.get(taxid)
        if rank in ncbi_rank_to_prefix:
            prefix = ncbi_rank_to_prefix[rank]
            name = names.get(taxid, "unknown")
            lineage_parts.append(f"{prefix}{name}")

    return "; ".join(lineage_parts)


def update_ncbi_taxonomy() -> NCBITaxa:
    return NCBITaxa()


def write_kraken_summary(
    report_counts: dict[str, dict[str, str]],
    fp: Path,
    ncbi: NCBITaxa,
) -> None:
    with open(fp, "w") as f:
        writer = csv.writer(f, delimiter="\t")
        header = ["#OTU ID"] + list(report_counts.keys()) + ["Consensus Lineage"]
        writer.writerow(header)

        all_ids = []
        for _, v in report_counts.items():
            all_ids.extend(v.keys())

        # Remove duplicates while preserving order
        all_ids = list(dict.fromkeys(all_ids))

        for taxon_id in all_ids:
            row = [taxon_id]
            for report_name in report_counts.keys():
                row.append(report_counts[report_name].get(taxon_id, "0"))
            row.append(get_consensus_lineage(taxon_id, ncbi))
            writer.writerow(row)


def main(reports: list[Path], summary: Path) -> None:
    report_counts = {
        report.stem.replace("-taxa", ""): parse_kraken(report) for report in reports
    }
    with tempfile.TemporaryDirectory() as tmpdir:
        fp = Path(tmpdir)
        os.environ["ETE_CACHE_DIR"] = fp.as_posix()
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
