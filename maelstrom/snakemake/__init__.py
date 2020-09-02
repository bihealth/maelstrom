"""maelstrom-snakemake -- call Snakemake with the Maelstrom Snakemake file."""

import pathlib
import shlex
import subprocess
import sys
import typing

#: Path to Snakefile.
PATH_SNAKEFILE = pathlib.Path(__file__).parent / "Snakefile"


def which(name: str) -> typing.Optional[str]:
    """Calls ``which`` and returns the path to the executable or None if there is none."""
    try:
        return (
            subprocess.check_output(["which", shlex.quote(name)], stderr=subprocess.PIPE)
            .decode("utf-8")
            .strip()
        )
    except subprocess.CalledProcessError:
        return None


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
