"""
Script for the download of files containing
the tables for the BDT trainings from hyperloop

\author Fabrizio Grosa <fabrizio.grosa@cern.ch>, CERN
\author Alexandre Bigot <alexandre.bigot@cern.ch>, Strasbourg University
\author Biao Zhang <biao.zhang@cern.ch>, CCNU
"""

import os
import argparse
from ROOT import TGrid  # pylint: disable=no-name-in-module


def main(infile, outpath):
    """
    Main function

    Parameters
    -----------------
    - infile: input file with list of directories separated by ','
    - outpath: output path
    """

    with open(infile) as f_txt:
        contents = f_txt.read()
    list_of_dirs = contents.split(",")

    grid = TGrid.Connect("alien://")

    for i_dir, indir in enumerate(list_of_dirs):
        jobdir = indir.split("/")[-1]
        outdir = os.path.join(outpath, jobdir)
        if not os.path.isdir(outdir):
            os.mkdir(outdir)
        list_of_subdirs = grid.Ls(indir.replace("alien://", ""))
        print(f"\033[32mDownload outputs from directory {i_dir+1}/{len(list_of_dirs)}"
              f" ({indir})\033[0m")
        is_merged = -1
        for i_sub, _ in enumerate(list_of_subdirs):
            if "AOD" in list_of_subdirs.GetFileName(i_sub):
                is_merged = i_sub
                break
        if is_merged >= 0:
            subdir = list_of_subdirs.GetFileName(is_merged)
            os.system("alien_cp -y 2 -T 32 -name contain_root "
                      f"{indir}/{subdir} file://{outdir}")
        else:
            for i_sub, _ in enumerate(list_of_subdirs):
                if list_of_subdirs.GetFileName(i_sub).isdigit():
                    subdir = list_of_subdirs.GetFileName(i_sub)
                    os.system("alien_cp -y 2 -T 32 -name contain_root "
                              f"{indir}/{subdir} file://{outdir}")


if __name__ == "__main__":
    PARSER = argparse.ArgumentParser(description="Arguments")
    PARSER.add_argument("infile", metavar="text", default="list.txt",
                        help="list of directories with input files")
    PARSER.add_argument("outpath", metavar="text", default=".",
                        help="output path")
    ARGS = PARSER.parse_args()

    main(ARGS.infile, ARGS.outpath)
