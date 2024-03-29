import os
import pytest
import shutil
import subprocess as sp
import sys
from pathlib import Path


@pytest.fixture
def setup(tmpdir):
    reads_fp = Path(".tests/data/reads/").resolve()
    hosts_fp = Path(".tests/data/hosts/").resolve()
    db_fp = Path(".tests/data/db/").resolve()

    project_dir = tmpdir / "project"

    sp.check_output(["sunbeam", "init", "--data_fp", reads_fp, project_dir])

    config_fp = project_dir / "sunbeam_config.yml"

    config_str = f"sbx_kraken: {{kraken_db_fp: {db_fp}}}"
    sp.check_output(
        [
            "sunbeam",
            "config",
            "modify",
            "-i",
            "-s",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    config_str = f"qc: {{host_fp: {hosts_fp}}}"
    sp.check_output(
        [
            "sunbeam",
            "config",
            "modify",
            "-i",
            "-s",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    yield tmpdir, project_dir


@pytest.fixture
def run_sunbeam(setup):
    tmpdir, project_dir = setup

    output_fp = project_dir / "sunbeam_output"

    try:
        # Run the test job
        sp.check_output(
            [
                "sunbeam",
                "run",
                "--conda-frontend",
                "conda",
                "--profile",
                project_dir,
                "all_classify",
                "--directory",
                tmpdir,
            ]
        )
    except sp.CalledProcessError as e:
        shutil.copytree(output_fp / "logs", "logs/")
        shutil.copytree(project_dir / "stats", "stats/")
        sys.exit(e)

    shutil.copytree(output_fp / "logs", "logs/")
    shutil.copytree(project_dir / "stats", "stats/")

    benchmarks_fp = project_dir / "stats"

    yield output_fp, benchmarks_fp


def test_full_run(run_sunbeam):
    output_fp, benchmarks_fp = run_sunbeam

    all_samples_fp = output_fp / "classify" / "kraken" / "all_samples.tsv"

    # Check output
    assert all_samples_fp.exists()

    with open(all_samples_fp) as f:
        header_line = f.readline()
        print(f"Header line: {header_line}")
        assert "TEST-taxa" in header_line
        assert "EMPTY-taxa" in header_line
        assert "Consensus Lineage" in header_line
        test_index = header_line.split("\t").index("TEST-taxa")
        empty_index = header_line.split("\t").index("EMPTY-taxa")

        lines = f.readlines()
        print(lines)
        for line in lines:
            if line[0] == "2":
                fields = line.split("\t")
                assert int(fields[empty_index]) == 0
                assert int(fields[test_index]) > 0
