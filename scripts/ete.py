import csv
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
