#!/usr/bin/env python3
import os
import logging
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from Bio import SeqIO

def dt_heatmap(concat_bed, bigwig, tag, mode, extra_args, output_dir, logger):
    """
    Helper for deepTools computeMatrix + plotHeatmap.
    Only called if bigwig is not None.
    """
    mat = os.path.join(output_dir, f"{tag}.mat.gz")
    img = os.path.join(output_dir, f"{tag}.png")
    if mode == 'rp':
        cmd = (
            f"computeMatrix reference-point --referencePoint center {extra_args} "
            f"-S {bigwig} -R {concat_bed} -o {mat} -p 80"
        )
    else:
        cmd = (
            f"computeMatrix scale-regions {extra_args} "
            f"-S {bigwig} -R {concat_bed} -o {mat} -p 80"
        )
    logger.info(cmd);  os.system(cmd)
    cmap_dict = {'insul':'Greens','chip':'Oranges_r','phast':'Blues_r','epi':'Greys'}
    cmap = cmap_dict.get(tag, 'Greys')
    plot_cmd = (
        f"plotHeatmap -m {mat} -out {img} "
        f"--colorMap {cmap} --whatToShow \"heatmap and colorbar\" "
        f"--sortRegions no --boxAroundHeatmaps no "
        f"--dpi 300 --heatmapHeight 100 --heatmapWidth 30"
    )
    logger.info(plot_cmd);  os.system(plot_cmd)
    logger.info(f"{tag} heatmap: {img}")

def main():
    p = argparse.ArgumentParser(
        description="Plot heatmaps for NPLB clusters and optional tracks"
    )
    p.add_argument("-n","--nplb_path",    required=True,
                   help="Path to nplb_clustered.bed")
    p.add_argument("-r","--reference_fasta", required=False, default=None,
                   help="Path to reference FASTA for getfasta + seq heatmaps")
    p.add_argument("-i","--insul_path",      required=False, default=None,
                   help="Path to insulation bigWig")
    p.add_argument("-c","--ctcf_chip_path",  required=False, default=None,
                   help="Path to CTCF ChIP-seq bigWig")
    p.add_argument("-p","--phastcons_path",  required=False, default=None,
                   help="Path to phastCons bigWig")
    p.add_argument("-e","--epigenomic_dir",  required=False, default=None,
                   help="Directory of other epigenomic .bw files")
    args = p.parse_args()

    # derive output_dir as subfolder of the NPLB BED file's directory
    nplb_path = args.nplb_path
    nplb_dir = os.path.dirname(nplb_path)
    output_dir = os.path.join(nplb_dir, "plot_heatmaps")
    os.makedirs(output_dir, exist_ok=True)

    # --- Logging ---
    log_file = os.path.join(output_dir, 'process.log')
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()]
    )
    logger = logging.getLogger()
    logger.info(f"Output directory: {output_dir}")

    # --- Read & split BED ---
    df = pd.read_csv(nplb_path, sep='\t', header=None)
    if 4 not in df.columns:
        logger.error("NPLB BED missing strand column"); return
    df[5] = df[4]  # duplicate for downstream
    clusters = sorted(df[3].unique())
    by_cluster = {c: df[df[3]==c].copy() for c in clusters}
    logger.info(f"Found {len(by_cluster)} clusters")

    # --- Per-cluster BEDs (with special flip+shift for cluster 11) ---
    beds = []
    for c in clusters:
        sub = by_cluster[c].sample(frac=1, random_state=42).copy()
        if c == 11:
            sub.loc[:,[4,5]] = sub.loc[:,[4,5]].replace({'+':'-','-':'+'})
            plus = sub[4]=='+'; minus=sub[4]=='-'
            sub.loc[plus,1]+=3; sub.loc[plus,2]+=3
            sub.loc[minus,1]-=3; sub.loc[minus,2]-=3
        path = os.path.join(output_dir, f"nplb_cluster_{c}.bed")
        sub.to_csv(path, sep='\t', header=False, index=False)
        beds.append(path)
        logger.info(f"Wrote cluster BED: {path}")

    # --- Concatenate ---
    concat = os.path.join(output_dir, 'nplb_clusters_concatenated.bed')
    with open(concat,'w') as out:
        for b in beds:
            out.write(open(b).read())
    logger.info(f"Concatenated BED: {concat}")

    # --- FASTA extraction & seq heatmaps (if requested) ---
    if args.reference_fasta:
        fasta_out = os.path.join(output_dir,'nplb_clusters_sequences.fa')
        cmd = (
            f"bedtools getfasta -fi {args.reference_fasta} "
            f"-bed {concat} -fo {fasta_out} -s"
        )
        logger.info(cmd); os.system(cmd)
        logger.info(f"FASTA→ {fasta_out}")

    # --- deepTools heatmaps for each provided track ---
    if args.insul_path:
        dt_heatmap(concat, args.insul_path, 'insul','rp','-a 50000 -b 50000',output_dir,logger)
    if args.ctcf_chip_path:
        dt_heatmap(concat, args.ctcf_chip_path,'chip','rp','-a 1000 -b 1000',output_dir,logger)
    if args.phastcons_path:
        dt_heatmap(concat, args.phastcons_path,'phast','scale','',output_dir,logger)

    if args.epigenomic_dir:
        for fn in sorted(os.listdir(args.epigenomic_dir)):
            if not fn.endswith('.bw'): continue
            bw = os.path.join(args.epigenomic_dir, fn)
            tag = os.path.splitext(fn)[0]
            dt_heatmap(concat,bw,tag,'rp','-a 1000 -b 1000',output_dir,logger)

    # --- Sequence heatmaps (only if FASTA was generated) ---
    if args.reference_fasta:
        records = list(SeqIO.parse(fasta_out,'fasta'))
        seqs = [str(r.seq) for r in records]
        data = []
        idx = 0
        for c in clusters:
            L = len(by_cluster[c])
            for _ in range(L):
                data.append((c, seqs[idx]))
                idx += 1
        df_seq = pd.DataFrame(data,columns=[0,2])

        # unclustered
        mapping = {'A':0,'C':1,'G':2,'T':3}
        mat = np.array([[mapping.get(nuc,np.nan) for nuc in s] for s in df_seq[2]])
        mask = np.ma.masked_invalid(mat)
        cmap = mcolors.ListedColormap(['green','blue','yellow','red'])
        cmap.set_bad('white')
        norm = mcolors.BoundaryNorm(np.arange(-0.5,4.5,1),4)
        plt.figure(figsize=(30,100))
        im = plt.imshow(mask,aspect='auto',cmap=cmap,norm=norm)
        plt.axis('off'); plt.colorbar(im,fraction=0.02,pad=0.02)
        plt.title('Unclustered Sequences',fontsize=40)
        plt.savefig(os.path.join(output_dir,'unclustered_heatmap.png'),dpi=300)
        plt.close()

        # upstream clustered (same ordering)
        plt.figure(figsize=(30,100))
        plt.imshow(mask,aspect='auto',cmap=cmap,norm=norm)
        plt.axis('off')
        plt.title('Motif Signal Upstream of CTCF Sites',fontsize=40)
        plt.savefig(os.path.join(output_dir,'upstream_clusters_heatmap.png'),dpi=300)
        plt.close()

    logger.info("Done.")

if __name__=="__main__":
    main()