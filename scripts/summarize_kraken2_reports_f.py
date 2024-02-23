from typing import Iterator, TextIO


allowed_ranks = {"D": "k", "P": "p", "C": "c", "O": "o", "F": "f", "G": "g", "S": "s"}
rank_order = {"D": 1, "P": 2, "C": 3, "O": 4, "F": 5, "G": 6, "S": 7}


def parse_kraken2_tsv_report(
    file_handler: TextIO,
) -> Iterator[dict[str, float | int | str]]:
    for line in file_handler:
        data = line.strip().split("\t")
        if len(data) == 6:
            (
                percentage,
                fragments_covered,
                fragments_assigned,
                rank,
                taxon_id,
                scientific_name,
            ) = data
            if rank in allowed_ranks.keys():
                yield {
                    "percentage": float(percentage),
                    "fragments_covered": int(fragments_covered),
                    "fragments_assigned": int(fragments_assigned),
                    "rank": rank,
                    "taxon_id": int(taxon_id),
                    "scientific_name": scientific_name,
                }


def consensus_lineage_str(rank_stack: list[str]) -> str:
    missing_ranks = [k for k, v in rank_order.items() if v > len(rank_stack)]
    rank_stack += [f"{allowed_ranks[r]}__" for r in missing_ranks]
    return "; ".join(rank_stack)


def create_kraken2_tsv_report(
    reports: list[Iterator[dict[str, float | int | str]]], report_names: list[str]
) -> tuple[dict[str, dict[int, int]], dict[int, str]]:
    consensus_lineages = {}
    report_counts = {}

    for report, report_name in zip(reports, report_names):
        rank_stack = []
        counts = {}

        for line in report:
            if line["rank"] in allowed_ranks.keys():
                # Update fragments assigned count
                counts[line["taxon_id"]] = line["fragments_assigned"]

                # Update rank stack
                if len(rank_stack) < rank_order[line["rank"]]:
                    rank_stack.append(
                        f"{allowed_ranks[line['rank']]}__{line['scientific_name'].lstrip()}"
                    )
                elif len(rank_stack) == rank_order[line["rank"]]:
                    rank_stack[-1] = (
                        f"{allowed_ranks[line['rank']]}__{line['scientific_name'].lstrip()}"
                    )
                else:
                    rank_stack = rank_stack[: rank_order[line["rank"]]]
                    rank_stack[-1] = (
                        f"{allowed_ranks[line['rank']]}__{line['scientific_name'].lstrip()}"
                    )

                # Update consensus lineages
                if line["taxon_id"] not in consensus_lineages:
                    consensus_lineages[line["taxon_id"]] = consensus_lineage_str(
                        rank_stack
                    )

        # Update report counts
        report_counts[report_name] = counts

    return report_counts, consensus_lineages


def write_kraken2_tsv_summary(
    report_counts: dict[str, dict[int, int]],
    consensus_lineages: dict[int, str],
    file_handler: TextIO,
) -> None:
    print(report_counts)
    print(consensus_lineages)
    # Write header
    file_handler.write(
        f"#OTU ID\t{'\t'.join([k for k, _ in report_counts.items()])}\tConsensus Lineage\n"
    )

    # Loop through consensus lineages
    for taxon_id, lineage in consensus_lineages.items():
        output_line = f"{taxon_id}\t"
        for report_name in [k for k, _ in report_counts.items()]:
            output_line += f"{report_counts[report_name].get(taxon_id, 0)}\t"
        file_handler.write(output_line + f"{lineage}\n")
