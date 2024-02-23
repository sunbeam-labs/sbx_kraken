import glob
import pytest
import sys
from pathlib import Path

from scripts.summarize_kraken2_reports_f import (
    parse_kraken2_tsv_report,
    create_kraken2_tsv_report,
    write_kraken2_tsv_summary,
)


@pytest.fixture
def reports():
    report_fp = Path(".tests/data/kraken2-outputs/").resolve()
    reports = glob.glob(str(report_fp / "*-taxa.tsv"))
    return reports, report_fp


def test_write_kraken2_tsv_summary(tmpdir, reports):
    reports, report_fp = reports
    report_names = [Path(x).stem for x in reports]

    report_counts, consensus_lineages = create_kraken2_tsv_report(
        [parse_kraken2_tsv_report(open(x)) for x in reports], report_names
    )
    print(report_counts)
    # print(consensus_lineages)

    summary_fp = tmpdir / "summary.tsv"
    summary_fp = report_fp / "summary.tsv"
    write_kraken2_tsv_summary(report_counts, consensus_lineages, open(summary_fp, "w"))

    with open(summary_fp, "r") as f:
        summary = f.readlines()
        report_generators = [parse_kraken2_tsv_report(open(x)) for x in reports]

        assert (
            summary[0].strip()
            == "#OTU ID\t" + "\t".join(report_names) + "\tConsensus Lineage"
        )
        assert (
            summary[1].strip()
            == "2\t"
            + "\t".join([str(next(x)["fragments_assigned"]) for x in report_generators])
            + "\tk__Bacteria; p__; c__; o__; f__; g__; s__"
        )


def test_parse_kraken2_tsv_report(reports):
    reports, _ = reports
    report = reports[0]
    parsed_report = parse_kraken2_tsv_report(open(report))

    assert list(next(parsed_report).keys()) == [
        "percentage",
        "fragments_covered",
        "fragments_assigned",
        "rank",
        "taxon_id",
        "scientific_name",
    ]
