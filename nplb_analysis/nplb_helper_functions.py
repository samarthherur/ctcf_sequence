import os
import argparse

def parse_nplb_args():
    parser = argparse.ArgumentParser(description="Run promoterLearn on neighbourhood FASTA")
    parser.add_argument("--fasta", "-f", required=True, help="Input FASTA file")
    parser.add_argument("--output_prefix", "-o", required=True, help="Prefix for output files")
    return parser.parse_args()

def run_promoter_learn(fasta_path: str, output_prefix: str):
    cmd = f"promoterLearn -f {fasta_path} -o {output_prefix}"
    print(f"Executing: {cmd}")
    ret = os.system(cmd)
    if ret != 0:
        raise RuntimeError(f"promoterLearn failed with exit code {ret}")