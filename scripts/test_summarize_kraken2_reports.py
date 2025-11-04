import os
import pytest
import sys
import tempfile
from pathlib import Path
from scripts.summarize_kraken2_reports import (
    allowed_ranks,
    get_consensus_lineage,
    parse_kraken,
    update_ncbi_taxonomy,
)


@pytest.fixture
def reports():
    report_fp = Path(".tests/data/kraken2-outputs/").resolve()
    return report_fp.glob("*-taxa.tsv")


@pytest.fixture
def test_update_consensus_lineage(autouse=True):
    with tempfile.TemporaryDirectory() as tmpdir:
        os.environ["ETE_CACHE_DIR"] = tmpdir
        ncbi = update_ncbi_taxonomy()


def test_parse_kraken(reports: list[Path]):
    for r in reports:
        k = parse_kraken(r)
        for i, j in k.items():
            int(i)  # taxon_id should be convertible to int
            int(j)  # fragments_assigned should be convertible to int


def test_get_consensus_lineage():
    ncbi = update_ncbi_taxonomy()
    lineage = get_consensus_lineage("562", ncbi)
    expected = "k__Bacteria; p__Pseudomonadota; c__Gammaproteobacteria; o__Enterobacterales; f__Enterobacteriaceae; g__Escherichia; s__Escherichia coli"
    assert lineage == expected, f"Expected: {expected}, Got: {lineage}"
    assert (
        get_consensus_lineage("000000", ncbi) == "unknown"
    ), "Unknown taxon ID should return 'unknown'"


def test_get_consensus_lineage_with_missing_family():
    ncbi = update_ncbi_taxonomy()
    lineage = get_consensus_lineage("754249", ncbi)
    expected = (
        "k__Bacteria; p__Actinomycetota; c__Actinomycetes; o__Micrococcales; "
        "f__; g__Luteimicrobium; s__"
    )
    assert lineage == expected, f"Expected: {expected}, Got: {lineage}"
