import os
import argparse
import pandas as pd

def parse_nplb_train_args():
    """
    Parse arguments for training: only input FASTA and output prefix.
    """
    parser = argparse.ArgumentParser(
        description="Run promoterLearn training on neighbourhood FASTA"
    )
    parser.add_argument(
        "--fasta", "-f",
        required=True,
        help="Input FASTA file (neighbourhood sequences) for training"
    )
    parser.add_argument(
        "--output_prefix", "-o",
        required=True,
        help="Prefix (directory+basename) for all training outputs"
    )
    return parser.parse_args()

def parse_nplb_classify_args():
    """
    Parse arguments for classification: input FASTA, model, output prefix,
    and a TSV mapping file with columns 'original_cluster_id', 'mapped_cluster_id', 'flip'.
    """
    parser = argparse.ArgumentParser(
        description="Run promoterClassify on neighbourhood FASTA with cluster mapping"
    )
    parser.add_argument(
        "--fasta", "-f",
        required=True,
        help="Input FASTA file (neighbourhood sequences) for classification"
    )
    parser.add_argument(
        "--model", "-m",
        required=True,
        help="Path to the trained promoterClassify model (.p file)"
    )
    parser.add_argument(
        "--output_prefix", "-o",
        required=True,
        help="Prefix (directory+basename) for all classification outputs"
    )
    parser.add_argument(
        "--cluster_map_tsv",
        required=True,
        help="TSV file mapping original_cluster_id to mapped_cluster_id, with a 'flip' column"
    )
    return parser.parse_args()

def run_promoter_learn(fasta_path: str, output_prefix: str):
    cmd = f"promoterLearn -f {fasta_path} -o {output_prefix}"
    print(f"Executing: {cmd}")
    ret = os.system(cmd)
    if ret != 0:
        raise RuntimeError(f"promoterLearn failed with exit code {ret}")

def run_promoter_classify(fasta_path: str, model_path: str, output_prefix: str):
    """
    Run promoterClassify on the given FASTA using a pre-trained model.
    """
    cmd = f"promoterClassify -f {fasta_path} -m {model_path} -o {output_prefix}"
    print(f"Executing: {cmd}")
    ret = os.system(cmd)
    if ret != 0:
        raise RuntimeError(f"promoterClassify failed with exit code {ret}")