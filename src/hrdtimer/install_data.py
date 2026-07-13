"""Post-install helper: download SigProfiler reference genomes.

Exposed as the ``hrdtimer-install-genomes`` console script. This is a one-time,
multi-GB download that cannot be performed by pip, so it is run explicitly after
installation.
"""

from __future__ import annotations

import argparse

from ._optional import require

DEFAULT_BUILDS = ("GRCh37", "GRCh38")


def install_genomes(builds=DEFAULT_BUILDS) -> None:
    """Install the given SigProfilerMatrixGenerator reference genome builds."""
    require("SigProfilerMatrixGenerator")  # clear error if the package is missing
    from SigProfilerMatrixGenerator import install as genInstall
    for build in builds:
        print(f" -- Installing reference genome {build} ...")
        genInstall.install(build, rsync=False, bash=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Install SigProfiler reference genomes for HRDTimer."
    )
    parser.add_argument(
        "builds",
        nargs="*",
        default=list(DEFAULT_BUILDS),
        help=f"Genome builds to install (default: {' '.join(DEFAULT_BUILDS)}).",
    )
    args = parser.parse_args()
    install_genomes(args.builds)


if __name__ == "__main__":
    main()
