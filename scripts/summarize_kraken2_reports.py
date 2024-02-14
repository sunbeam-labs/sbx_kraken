from collections import OrderedDict
from pathlib import Path
from typing import Iterator, TextIO

def parse_kraken2_tsv_report(file_handler: TextIO) -> Iterator[dict[str, float | int | str]]:
    for line in file_handler:
        data = line.strip().split('\t')
        if len(data) == 6:
            percentage, fragments_covered, fragments_assigned, rank, taxon_id, scientific_name = data
            yield {
                'percentage': float(percentage),
                'fragments_covered': int(fragments_covered),
                'fragments_assigned': int(fragments_assigned),
                'rank': rank,
                'taxon_id': int(taxon_id),
                'scientific_name': scientific_name
            }

def write_kraken2_tsv_summary(reports: list[Iterator[dict[str, float | int | str]]], report_names: list[str], file_handler: TextIO) -> None:
    next_lines = OrderedDict(zip(report_names, [next(r, {}) for r in reports]))

    # Write header
    file_handler.write(f"#OTU ID\t{'\t'.join(report_names)}\n")

    # Loop through report generators until all are exhausted
    while any([r for _, r in next_lines.items()]):
        # Check indentation levels of each scientific_name to determine rank
        # Collect list of lowest rank lines
        # Only keep the one(s) with the lowest OTU ID
        lowest_rank_lines = {}
        lowest_rank = float('inf')
        for n, r in next_lines.items():
            if r:
                if len(r['scientific_name']) - len(r['scientific_name'].lstrip()) < lowest_rank:
                    lowest_rank = len(r['scientific_name']) - len(r['scientific_name'].lstrip(' '))
                    lowest_rank_lines = {n: r}
                elif len(r['scientific_name']) - len(r['scientific_name'].lstrip()) == lowest_rank:
                    lowest_rank_lines[n] = r

        # Check that the lowest rank lines have the same OTU ID
        otu_set = set([r['taxon_id'] for r in lowest_rank_lines.values()])
        assert len(otu_set) == 1, "OTU IDs do not match for lowest rank lines"
        otu_id = list(otu_set)[0]

        # Write the line
        output_line = f"{otu_id}\t"
        for n, r in next_lines.items():
            if r['taxon_id'] == otu_id:
                output_line += f"{r['fragments_assigned']}\t"
                next_lines[n] = next(reports[report_names.index(n)], {})
            else:
                output_line += "0\t"



report_names = [Path(x).stem for x in snakemake.input.reports]
write_kraken2_tsv_summary([parse_kraken2_tsv_report(open(x)) for x in snakemake.input.reports], report_names, open(snakemake.output.summary))