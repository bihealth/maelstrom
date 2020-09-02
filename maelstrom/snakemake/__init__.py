"""maelstrom-snakemake -- call Snakemake with the Maelstrom Snakemake file."""

import pathlib
import subprocess
import sys
import typing

#: Path to Snakefile.
PATH_SNAKEFILE = pathlib.Path(__file__).parent / "Snakefile"


def main(args=None) -> int:
    try:
        subprocess.check_call(
            ["snakemake", "--snakefile", str(PATH_SNAKEFILE), *(args or sys.argv[1:])]
        )
    except subprocess.CalledProcessError:
        return 1
    else:
        return 0


if __name__ == "__main__":
    sys.exit(main())
