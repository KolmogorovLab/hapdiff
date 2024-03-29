#!/usr/bin/env python3

from threading import Thread
import sys
import subprocess
import os
import argparse
from distutils import spawn

from Bio import SeqIO
from Bio.SeqIO import SeqRecord


pipeline_dir = os.path.dirname(os.path.realpath(__file__))
MINIMAP2 = os.path.join(pipeline_dir, "submodules", "minimap2", "minimap2")
SAMTOOLS = "samtools"
BEDTOOLS = "bedtools"

VERSION = "0.9"

sys.path.insert(0, os.path.join(pipeline_dir, "submodules", "svim-asm", "src"))
import svim_asm.main as svim


def file_check(filename):
    if not os.path.isfile(filename) or os.path.getsize(filename) == 0:
        raise Exception("File not found, or has zero length:", filename)


def generate_alignment(ref_path, asm_path, num_threads, out_bam):
    cmd = MINIMAP2 + " -ax asm20 -B 2 -E 3,1 -O 6,100 --cs -t {0} {1} {2} -K 5G | samtools sort -m 4G -@ 8 >{3}" \
                            .format(num_threads, ref_path, asm_path, out_bam)
    print("Running: " + cmd)
    subprocess.check_call(cmd, shell=True)
    subprocess.check_call("samtools index -@ 4 {0}".format(out_bam), shell=True)


def fragment(input_fasta, output_fasta, frag_size):
    with open(output_fasta, "w") as fout:
        for seq in SeqIO.parse(input_fasta, "fasta"):
            if len(seq.seq) < frag_size:
                SeqIO.write(seq, fout, "fasta")
            else:
                for chunk in range(0, len(seq.seq) // frag_size + 1):
                    chunk_seq = seq.seq[chunk * frag_size : (chunk + 1) * frag_size]
                    chunk_id = str(seq.id) + "_chunk_" + str(chunk)
                    SeqIO.write(SeqRecord(seq=chunk_seq, id=chunk_id, description=""), fout, "fasta")


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
    parser.add_argument("--tandem-repeats", dest="tandem_repeats",
                        default=None, required=False,
                        metavar="path", help="Tandem repeat annotation in bed format")
    parser.add_argument("--sample", dest="sample",
                        default="Sample", required=False,
                        help="Sample ID [deafult=Fample]")
    parser.add_argument("--sv-size", dest="sv_size", type=int,
                        default=30, metavar="int", help="minimum SV size [30]")
    parser.add_argument("--fragment", dest="fragment", type=int,
                        default=None, metavar="int", help="fragment query to X Mb to reduce minimap2 memory footprint [None]")
    #parser.add_argument("--phased", dest="phased", action="store_true",
    #                    default=False, help="produce phased vcf")
    parser.add_argument("-t", "--threads", dest="threads", type=int,
                        default=10, metavar="int", help="number of parallel threads [10]")
    parser.add_argument("-v", "--version", action="version", version=VERSION)
    args = parser.parse_args()

    for e in [MINIMAP2, SAMTOOLS, BEDTOOLS]:
        if not spawn.find_executable(e):
            print("Not installed: " + e, file=sys.stderr)
            return 1

    if not os.path.isdir(args.out_dir):
        os.mkdir(args.out_dir)

    file_check(args.reference)
    file_check(args.hap_pat)
    file_check(args.hap_mat)

    prefix = "hapdiff"
    aln_1 = os.path.join(args.out_dir, prefix + "_pat" + ".bam")
    aln_2 = os.path.join(args.out_dir, prefix + "_mat" + ".bam")

    fragmented_pat = args.hap_pat
    fragmented_mat = args.hap_mat
    if args.fragment is not None:
        fragmented_pat = os.path.join(args.out_dir, "fragmented_pat.fasta")
        fragment(args.hap_pat, fragmented_pat, args.fragment * 1000000)
        fragmented_mat = os.path.join(args.out_dir, "fragmented_mat.fasta")
        fragment(args.hap_mat, fragmented_mat, args.fragment * 1000000)

    generate_alignment(args.reference, fragmented_pat, args.threads, aln_1)
    file_check(aln_1)
    generate_alignment(args.reference, fragmented_mat, args.threads, aln_2)
    file_check(aln_2)

    def run_svim(out_file, phased):
        svim_cmd = ["diploid", args.out_dir, aln_1, aln_2, args.reference, "--min_sv_size", str(args.sv_size),
                    "--partition_max_distance", "5000", "--max_edit_distance", "0.3",
                    "--filter_contained", "--query_names", "--sample", args.sample]
        if phased:
            svim_cmd.append("--phased_gt")
        if args.tandem_repeats:
            svim_cmd.extend(["--tandem", args.tandem_repeats])
        svim.main(svim_cmd)

        SVIM_OUTPUT = os.path.join(args.out_dir, "variants.vcf")
        SV_LENGTHS = os.path.join(args.out_dir, "sv-lengths.png")
        out_with_prefix = os.path.join(args.out_dir, out_file)
        os.rename(SVIM_OUTPUT, out_with_prefix)

        subprocess.check_call(["bgzip", "-f", out_with_prefix])
        subprocess.check_call(["tabix", "-f", out_with_prefix + ".gz"])
        if os.path.isfile(SV_LENGTHS):
            os.remove(SV_LENGTHS)

    run_svim("hapdiff_unphased.vcf", False)
    run_svim("hapdiff_phased.vcf", True)

    conf_pat = os.path.join(args.out_dir, "aln_coverage_pat.bed")
    conf_mat = os.path.join(args.out_dir, "aln_coverage_mat.bed")
    merged_bed = os.path.join(args.out_dir, "confident_regions.bed")

    with open(merged_bed, "w") as fout:
        fout.write("#SAMPLE:{0}\n".format(args.sample))
        fout.flush()
        bedtools_cmd = [BEDTOOLS, "intersect", "-a", conf_pat, "-b", conf_mat, "-sortout", "|", "uniq"]
        subprocess.check_call(" ".join(bedtools_cmd), shell=True, stdout=fout)

    return 0


if __name__ == "__main__":
    sys.exit(main())
