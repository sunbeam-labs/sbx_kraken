import os
import shutil
import subprocess as sp
import tempfile
import unittest


class FullRunTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp_dir = tempfile.mkdtemp()

        self.reads_fp = os.path.abspath(".tests/data/reads/")
        self.db_fp = os.path.abspath(".tests/data/db/")

        self.project_dir = os.path.join(self.temp_dir, "project/")

        sp.check_output(
            ["sunbeam", "init", "--data_fp", self.reads_fp, self.project_dir]
        )

        self.config_fp = os.path.join(self.project_dir, "sunbeam_config.yml")

        config_str = f"sbx_kraken: {{kraken_db_fp: {self.db_fp}}}"

        sp.check_output(
            [
                "sunbeam",
                "config",
                "modify",
                "-i",
                "-s",
                f"{config_str}",
                f"{self.config_fp}",
            ]
        )

        self.output_fp = os.path.join(self.project_dir, "sunbeam_output")
        # shutil.copytree(".tests/data/sunbeam_output", self.output_fp)

        self.all_samples_fp = os.path.join(
            self.output_fp, "classify/kraken/all_samples.tsv"
        )

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_full_run(self):
        # Run the test job.
        sp.check_output(
            [
                "sunbeam",
                "run",
                "--profile",
                self.project_dir,
                "all_classify",
                "--directory",
                self.temp_dir,
            ]
        )

        # Check output
        self.assertTrue(os.path.exists(self.all_samples_fp))
        with open(self.all_ptr_fp) as f:
            self.assertEqual(next(f), "\tTEST0\tTEST1\tTEST2\tTEST3\tTEST4\n")
            for val in next(f).split("\t")[1:]:
                self.assertEqual(round(float(val)), 3)
