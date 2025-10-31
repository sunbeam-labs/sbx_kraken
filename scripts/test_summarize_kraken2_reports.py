import glob
import pytest
import sys
from pathlib import Path

from .summarize_kraken2_reports import (
    parse_kraken,
    get_consensus_lineage,
    write_kraken_summary,
)


@pytest.fixture
def reports():
    report_fp = Path(".tests/data/kraken2-outputs/").resolve()
    reports = glob.glob(str(report_fp / "*-taxa.tsv"))
    return reports, report_fp