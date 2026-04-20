from pathlib import Path
from typing import Iterator, TextIO


def parse_kraken2_tsv_report(file_handler: TextIO) -> Iterator[dict]:
    for line in file_handler:
        parts = line.rstrip("\n").split("\t")
        if len(parts) != 6:
            continue
        percent, clade_reads, direct_reads, rank_code, taxid, name = parts
        depth = len(name) - len(name.lstrip())  # leading spaces as taxonomic depth
        yield {
            "percentage": float(percent),
            "fragments_covered": int(clade_reads),
            "fragments_assigned": int(direct_reads),
            "rank_code": rank_code,
            "taxon_id": int(taxid),
            "scientific_name": name.strip(),
            "depth": depth // 2,  # kraken2 uses 2 spaces per level
        }


def consensus_lineage_str(lineage_stack: list[str]) -> str:
    # Fill missing levels with empty placeholders if needed
    full_lineage = lineage_stack + [f"__"] * (max(7 - len(lineage_stack), 0))
    return "; ".join(full_lineage)


def create_kraken2_tsv_report(
    reports: list[Iterator[dict]], report_names: list[str]
) -> tuple[dict[str, dict[int, int]], dict[int, str]]:
    consensus_lineages = {}
    report_counts = {}

    for report, report_name in zip(reports, report_names):
        lineage_stack = []
        counts = {}

        for entry in report:
            depth = entry["depth"]
            name = entry["scientific_name"]
            taxon_id = entry["taxon_id"]
            assigned = entry["fragments_assigned"]

            # Update stack for lineage at current depth
            if len(lineage_stack) <= depth:
                lineage_stack.extend(["__"] * (depth - len(lineage_stack) + 1))
            lineage_stack = lineage_stack[: depth + 1]
            lineage_stack[depth] = f"__{name}"

            # Store read count
            counts[taxon_id] = assigned

            # Store consensus lineage string (once)
            if taxon_id not in consensus_lineages:
                consensus_lineages[taxon_id] = consensus_lineage_str(lineage_stack)

        report_counts[report_name] = counts

    return report_counts, consensus_lineages


def write_kraken2_tsv_summary(
    report_counts: dict[str, dict[int, int]],
    consensus_lineages: dict[int, str],
    file_handler: TextIO,
) -> None:
    sample_names = list(report_counts.keys())
    header = "#OTU ID\t" + "\t".join(sample_names) + "\tConsensus Lineage\n"
    file_handler.write(header)

    for taxid, lineage in consensus_lineages.items():
        counts = [str(report_counts[sample].get(taxid, 0)) for sample in sample_names]
        file_handler.write(f"{taxid}\t" + "\t".join(counts) + f"\t{lineage}\n")


# If using Snakemake
if "snakemake" in globals():
    reports = snakemake.input.reports  # type: ignore
    summary = snakemake.output.summary  # type: ignore

    report_names = [Path(p).stem for p in reports]
    parsed_reports = [parse_kraken2_tsv_report(open(p)) for p in reports]

    report_counts, consensus_lineages = create_kraken2_tsv_report(
        parsed_reports, report_names
    )

    with open(summary, "w") as out_f:
        write_kraken2_tsv_summary(report_counts, consensus_lineages, out_f)
