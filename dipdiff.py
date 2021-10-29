#!/usr/bin/env python3

from threading import Thread
import sys
import subprocess
import os
import argparse
from distutils import spawn


pipeline_dir = os.path.dirname(os.path.realpath(__file__))
MINIMAP2 = os.path.join(pipeline_dir, "submodules", "minimap2", "minimap2")
SAMTOOLS = "samtools"

sys.path.insert(0, os.path.join(pipeline_dir, "submodules", "svim-asm", "src"))
import svim_asm.main as svim


def file_check(filename):
    if not os.path.isfile(filename) or os.path.getsize(filename) == 0:
        raise Exception("File not found, or has zero length:", filename)


def generate_alignment(ref_path, asm_path, num_threads, out_bam):
    cmd = "minimap2 -ax asm20 -B 2 -E 3,1 -O 6,100 --cs -t {0} {1} {2} -K 5G | samtools sort -m 4G -@ 8 >{3}" \
                            .format(num_threads, ref_path, asm_path, out_bam)
    print("Running: " + cmd)
    subprocess.check_call(cmd, shell=True, stderr=open(os.devnull, "w"))
    subprocess.check_call("samtools index -@ 4 {0}".format(out_bam), shell=True)


def main():
    if sys.version_info < (3,):
        raise SystemExit("Requires Python 3")

    parser = argparse.ArgumentParser \
        (description="Call structural variants for a diploid assembly")

    parser.add_argument("--reference", dest="reference",
                        metavar="path", required=True,
                        help="path to reference file (fasta format)")
    parser.add_argument("--pat", dest="hap_pat", required=True, metavar="path",
                        help="path to paternal haplotype (in fasta format)")
    parser.add_argument("--mat", dest="hap_mat", required=True, metavar="path",
                        help="path to maternal haplotype (in fasta format)")
    parser.add_argument("--out-dir", dest="out_dir",
                        default=None, required=True,
                        metavar="path", help="Output directory")
    parser.add_argument("--sv-size", dest="sv_size", type=int,
                        default=50, metavar="int", help="minimum SV size [50]")
    parser.add_argument("-t", "--threads", dest="threads", type=int,
                        default=10, metavar="int", help="number of parallel threads [10]")
    args = parser.parse_args()

    for e in [MINIMAP2, SAMTOOLS]:
        if not spawn.find_executable(e):
            print("Not installed: " + e, file=sys.stderr)
            return 1

    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)

    file_check(args.reference)
    file_check(args.hap_pat)
    file_check(args.hap_mat)

    prefix = "dipdiff"
    aln_1 = os.path.join(args.out_dir, prefix + "_pat" + ".bam")
    aln_2 = os.path.join(args.out_dir, prefix + "_mat" + ".bam")

    thread_1 = Thread(target=generate_alignment, args=(args.reference, args.hap_pat, args.threads // 2, aln_1))
    thread_2 = Thread(target=generate_alignment, args=(args.reference, args.hap_mat, args.threads // 2, aln_2))
    thread_1.start()
    thread_2.start()
    thread_1.join()
    thread_2.join()

    file_check(aln_1)
    file_check(aln_2)

    svim_cmd = ["diploid", args.out_dir, aln_1, aln_2, args.reference, "--min_sv_size", str(args.sv_size),
                "--partition_max_distance", "5000", "--max_edit_distance", "0.3", "--filter_contained"]
    svim.main(svim_cmd)
    #os.remove(os.path.join(out_dir, "sv-lengths.png"))

    return 0


if __name__ == "__main__":
    sys.exit(main())
