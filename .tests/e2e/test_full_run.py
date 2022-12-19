import csv
import os
import pytest
import shutil
import subprocess as sp
import tempfile


@pytest.fixture
def setup():
    temp_dir = tempfile.mkdtemp()

    reads_fp = os.path.abspath(".tests/data/reads/")
    db_fp = os.path.abspath(".tests/data/db/")

    project_dir = os.path.join(temp_dir, "project/")

    sp.check_output(["sunbeam", "init", "--data_fp", reads_fp, project_dir])

    config_fp = os.path.join(project_dir, "sunbeam_config.yml")

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

    yield temp_dir, project_dir

    shutil.rmtree(temp_dir)


@pytest.fixture
def run_sunbeam(setup):
    temp_dir, project_dir = setup

    # Run the test job.
    sp.check_output(
        [
            "sunbeam",
            "run",
            "--profile",
            project_dir,
            "all_classify",
            "--directory",
            temp_dir,
        ]
    )

    output_fp = os.path.join(project_dir, "sunbeam_output")

    all_samples_fp = os.path.join(output_fp, "classify/kraken/all_samples.tsv")

    benchmarks_fp = os.path.join(project_dir, "stats/")

    yield all_samples_fp, benchmarks_fp


def test_full_run(run_sunbeam):
    all_samples_fp, benchmarks_fp = run_sunbeam

    # Check output
    assert os.path.exists(all_samples_fp)


def test_benchmarks(run_sunbeam):
    all_samples_fp, benchmarks_fp = run_sunbeam

    filename = os.listdir(benchmarks_fp)[0]
    with open(os.path.join(benchmarks_fp, filename)) as f:
        rd = csv.DictReader(f, delimiter="\t")
        for r in rd:
            if r["rule"] in ["classic_k2_biom", "kraken2_biom"]:
                assert (
                    float(r["cpu_time"]) < 1
                ), f"cpu_time for {r['rule']} is higher than 1: {r['cpu_time']}"
            else:
                assert (
                    float(r["cpu_time"]) < 0.5
                ), f"cpu_time for {r['rule']} is higher than 0.5: {r['cpu_time']}"
