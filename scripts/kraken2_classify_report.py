import gzip
import subprocess as sp
from pathlib import Path

def not_empty():
    with gzip.open(snakemake.input[0], "rb") as f:
        try:
            file_content = f.read(1)
            return len(file_content) > 0
        except:
            return False


if not_empty():
    args = [
        "kraken2",
        "--gzip-compressed",
        "--db",
        f"{snakemake.params.db}",
        "--report",
        f"{snakemake.output.report}",
        "--output",
        f"{snakemake.output.raw}",
        f"{snakemake.params.paired_end}",
    ]
    [args.append(x) for x in snakemake.input]
    output = sp.check_output(args)

    with open(snakemake.log[0], "wb") as f:
        f.write(output)
else:
    with open(snakemake.log[0], "w") as f_log, open(
        snakemake.output.report, "w"
    ) as f_out:
        f_log.write("Empty reads files")
        f_out.write("0\t0.0\tk__Bacteria; p__; c__; o__; f__; g__; s__")
    Path(snakemake.output.raw).touch()