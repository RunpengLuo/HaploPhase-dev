import os
import sys
from lib.utils import *
from lib.file_io import *

def get_ctg_coverage(bam_file: str, out_dir: str):
    out_file = out_dir + "/contigs.coverage.txt"
    if not os.path.exists(bam_file + ".bai"):
        # generate the index file first
        System(f"samtools index {bam_file}")
    System(f"samtools coverage {bam_file} > {out_file}")

    return

def load_coverage_file(cov_file: str):
    cov_dict = {}
    with open(cov_file, "r") as cov_fd:
        cov_fd.readline()
        for entry in cov_fd:
            rname,_,_,_,_,coverage,meandepth,_,meanmapq = entry.strip().split()
            cov_dict[rname] = [float(coverage), float(meandepth), float(meanmapq)]
        cov_fd.close()
    return cov_dict