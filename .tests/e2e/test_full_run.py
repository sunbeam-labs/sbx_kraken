import os
import pytest
import shutil
import subprocess as sp
from pathlib import Path


@pytest.fixture
def setup(tmp_path):
    reads_fp = Path(".tests/data/reads/").resolve()
    db_fp = Path(".tests/data/db/").resolve()
    hosts_fp = Path(".tests/data/hosts/").resolve()

    project_dir = tmp_path / "project/"

    sp.check_output(["sunbeam", "init", "--data_fp", reads_fp, project_dir])

    config_fp = project_dir / "sunbeam_config.yml"

    config_str = f"sbx_kraken: {{kraken_db_fp: {db_fp}}}"
    sp.check_output(
        [
            "sunbeam",
            "config",
            "--modify",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    config_str = f"qc: {{host_fp: {hosts_fp}}}"
    sp.check_output(
        [
            "sunbeam",
            "config",
            "--modify",
            f"{config_str}",
            f"{config_fp}",
        ]
    )

    yield tmp_path, project_dir

    shutil.rmtree(tmp_path)


@pytest.fixture
def run_sunbeam(setup):
    tmp_path, project_dir = setup
    output_fp = project_dir / "sunbeam_output"
    log_fp = output_fp / "logs"
    stats_fp = project_dir / "stats"

    sbx_proc = sp.run(
        [
            "sunbeam",
            "run",
            "--profile",
            project_dir,
            "all_kraken",
            "--directory",
            tmp_path,
        ],
        capture_output=True,
        text=True,
    )

    print("STDOUT: ", sbx_proc.stdout)
    print("STDERR: ", sbx_proc.stderr)

    if os.getenv("GITHUB_ACTIONS") == "true":
        try:
            shutil.copytree(log_fp, "logs/")
            shutil.copytree(stats_fp, "stats/")
        except FileNotFoundError:
            print("No logs or stats directory found.")

    output_fp = project_dir / "sunbeam_output"
    benchmarks_fp = project_dir / "stats/"

    yield output_fp, benchmarks_fp, sbx_proc


def test_full_run(run_sunbeam):
    output_fp, benchmarks_fp, proc = run_sunbeam

    assert proc.returncode == 0, f"Sunbeam run failed with error: {proc.stderr}"

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
