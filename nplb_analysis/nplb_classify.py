import os
from nplb_helper_functions import parse_nplb_args, run_promoter_learn

def main():
    args = parse_nplb_args()

    # 1) Run promoterLearn
    run_promoter_learn(args.fasta, args.output_prefix)

if __name__ == "__main__":
    main()