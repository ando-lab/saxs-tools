"""
Print the version number
"""

import argparse
import os
import sys

import saxs_tools

parser = argparse.ArgumentParser(
    description=__doc__,
)


def run(args=None):
    args = parser.parse_args(args)
    print("saxs-tools:", saxs_tools.__version__)
    print("Python {0.major}.{0.minor}.{0.micro}".format(sys.version_info))
    print(f"Installed in: {os.path.split(saxs_tools.__file__)[0]}")


if __name__ == "__main__":
    run()
