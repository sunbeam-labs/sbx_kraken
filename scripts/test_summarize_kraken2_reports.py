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

    summary_fp = tmpdir / "summary.tsv"
    summary_fp = report_fp / "summary.tsv"
    write_kraken2_tsv_summary(report_counts, consensus_lineages, open(summary_fp, "w"))

    with open(summary_fp, "r") as f:
        summary = f.readlines()
        summary = {x.split("\t")[0]: x.split("\t")[1:] for x in summary}
        report_generators = [parse_kraken2_tsv_report(open(x)) for x in reports]

        assert summary["#OTU ID"] == report_names + ["Consensus Lineage\n"]
        assert summary["2"] == [
            str(next(x)["fragments_assigned"]) for x in report_generators
        ] + ["k__Bacteria; p__; c__; o__; f__; g__; s__\n"]

        for rg in report_generators:
            for line in rg:
                assert str(line["fragments_assigned"]) in summary[str(line["taxon_id"])]
                assert (
                    str(line["scientific_name"]).lstrip()
                    in summary[str(line["taxon_id"])][-1]
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
