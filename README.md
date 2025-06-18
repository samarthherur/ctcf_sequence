## 1. Conceptual & Usage Overview of `~/ctcf_sequence/fimo_analysis/fimo_neighbourhood_analysis.py`

This script implements an end-to-end pipeline to take raw FIMO motif calls and integrate them with ChIP-seq and boundary data, producing statistical summaries and visualizations at each step.

1. **Inputs**  
   - **FIMO output** (TSV of motif hits)  
   - **Genome FASTA** (reference sequence)  
   - **ChIP-seq peaks** (BED file)  
   - **Boundary regions** (BED file)  
   - **Neighborhood size** (window around each motif, e.g. 50 bp)  
   - **Output directory** (where all results are saved)

2. **Neighborhood Extraction**  
   1. **Extend** each FIMO hit by ±_neighborhood_size_ around its genomic coordinates.  
   2. **Merge** overlapping windows so each motif has a unique “neighborhood” window.  
   3. **Split** by strand and write plus/minus BED files.  
   4. **Fetch** the corresponding sequences with `bedtools getfasta`.  
   5. **Save** a combined FASTA of all neighborhoods.

3. **Visualization & Initial Statistics**  
   - **Heatmap & Logo** of the raw neighborhood FASTA.  
   - **Median score** of all neighborhood windows.

4. **ChIP-seq Filtering**  
   - **Intersect** neighborhood windows with ChIP-seq peaks.  
   - **Save** filtered BED & FASTA.  
   - **Redo** heatmap & logo on filtered sequences.

5. **Boundary Intersection**  
   - **Split** filtered windows into:  
     - Those that **intersect** boundary regions  
     - Those that **do not** intersect  
   - **Heatmap & logo** for both subsets.

6. **Score Distributions & Tests**  
   - **Compute** median scores for intersecting vs. non-intersecting windows.  
   - **Plot** and save:  
     - Density overlay (`score_distribution.png`)  
     - Overlaid histograms (`score_histogram.png`)  
   - **Kolmogorov–Smirnov test** between the two groups.

7. **Motif Classification**  
   - **Filter** the original FIMO TSV by ChIP peaks.  
   - **Assign** each motif to boundary windows.  
   - **Classify** motifs into:  
     1. Class 1 (outside any boundary)  
     2. Class 2 (inside a boundary, no special neighbor relationship)  
     3. Class 3 (strand-paired motifs across adjacent boundaries)  
   - **Write** `classified_motifs.tsv`.

8. **Classified Motif Summaries**  
   - **Total count** and **per-class counts**.  
   - **Median score** per class.  
   - **KS tests** between every pair of classes.  
   - **Density plot** of class-specific score distributions.

9. **Repeat-Filtering & Final Visuals**  
   - **Remove** any sequences with `N`/`n` or lowercase letters.  
   - **Report** number of repeat-free sequences and per-base counts.  
   - **Save** a `_wr.fasta` of cleaned sequences.  
   - **Heatmap & logo** of repeat-free FASTA.

10. **Neighborhood Matrix**  
    - **Count** overlaps with ChIP and boundary for each window (`-c`).  
    - **Save** `neighbourhood_matrix.tsv` for downstream analyses.

---

**Usage Overview: Input Arguments**

| Argument               | Type       | Required  | Description                                                                                                                               |
|------------------------|------------|-----------|-------------------------------------------------------------------------------------------------------------------------------------------|
| `--fimo_filepath`      | string     | yes       | Path to the **raw FIMO TSV** file containing motif occurrences (chromosome, start, end, score, p-value, etc.).                            |
| `--genome_fasta`       | string     | yes       | Path to the **reference genome FASTA**. Used by `bedtools getfasta` to extract the ±neighborhood window around each motif hit.           |
| `--chip_filepath`      | string     | yes       | Path to the **ChIP-seq peaks BED** file. Windows that overlap these regions are retained during the ChIP-filtering step.                   |
| `--boundary`           | string     | yes       | Path to the **boundary regions BED** file. Filtered windows are split into “intersect” vs. “non-intersect” groups based on these coords.    |
| `--output_dir`         | string     | yes       | Directory where all outputs will be written (created if it doesn’t exist). Contains FASTA/BED files, PNG figures, matrices, and reports.  |
| `--neighbourhood_size` | integer    | no (50)   | Number of base-pairs to extend **upstream and downstream** of each motif hit when constructing its “neighborhood” window. Defaults to 50 bp. |
| `--ignore_repeats`     | flag       | no        | If provided, **skip** the step that filters out sequences containing `N`/`n` or lowercase letters from the neighborhood FASTA.            |

**Example Invocation**

```bash
python fimo_neighbourhood_analysis.py \
  --fimo_filepath   ~/data/fimo_ctcf.tsv \
  --genome_fasta    ~/genomes/hg38.fa \
  --chip_filepath   ~/data/ENCFF294RSZ.bed \
  --boundary        ~/data/boundaries.bed \
  --output_dir      ~/results/fimo_analysis/ \
  --neighbourhood_size 100 \
  --ignore_repeats
```

## 2. Conceptual & Usage Overview of `train.sh`

This Bash script is a SLURM batch submission that runs the NPLB `promoterLearn` training on a pre-filtered FIMO neighbourhood FASTA, configuring compute resources, environment modules, and logging.

### SLURM Directives  
- `#SBATCH --job-name=nplb_train`  
- `#SBATCH --output=nplb_train.out`  
- `#SBATCH --error=nplb_train.err`  
- `#SBATCH --time=5-10:00:00`  
- `#SBATCH --partition=gpu`  
- `#SBATCH --ntasks=40`  
- `#SBATCH --nodelist=cn1`  

---

### Inputs

| Flag | Parameter | Type   | Required | Description                                                      |
|------|-----------|--------|----------|------------------------------------------------------------------|
| `-f`  | `<fasta>` | string | yes      | Path to the filtered, repeat-free neighbourhood FASTA file.      |
| `-o`  | `<out_dir>` | string | yes      | Output directory for NPLB results (models, clustering outputs).  |

---

### Usage

Submit the job to SLURM:

```bash
sbatch train.sh
```
## 3. Conceptual & Usage Overview of `nplb_train.py`

This Python script parses the NPLB `promoterLearn` training output (`architectureDetails.txt`) and converts it into a BED file of clustered windows.

### Inputs

| Flag                           | Parameter           | Type   | Required | Description                                                                                               |
|--------------------------------|---------------------|--------|----------|-----------------------------------------------------------------------------------------------------------|
| `--nplb_train_dir`, `-o`        | `<nplb_train_dir>`   | string | yes      | Directory of NPLB training results containing `architectureDetails.txt`.     |

### Functions Utilized

- `parse_nplb_train_args()`: parses the `--nplb_train_dir` argument.  
  *(defined in `nplb_helper_functions.py`)*

### Processing Steps

1. **Locate** `architectureDetails.txt` in the directory of `nplb_train_dir`.  
2. **Assert** that the file exists, otherwise exit with an error message.  
3. **Load** the TSV into a pandas DataFrame (`nplb_clustered`).  
4. **Extract** from each entry string:  
   - Chromosome (`chr`)  
   - Start coordinate (`start`)  
   - End coordinate (`end`)  
   - Strand (`+` or `-`)  
   - Cluster ID (`cluster`)  
5. **Assemble** a new DataFrame with columns `['chr','start','end','cluster','strand']`.  
6. **Write** this DataFrame as a BED file named `nplb_clustered.bed` alongside the input file.  
7. **Print** the path of the saved BED.

### Usage

```bash
python nplb_train.py \
  --nplb_train_dir /home/samarth/ctcf_sequence_data/output/NPLB/HFF
```

## 4. Conceptual & Usage Overview of `nplb_ordering.py`

This Python script performs a multi‐step analysis on NPLB‐clustered windows, optionally computing  
- **Gene‐TSS distances** and hypergeometric enrichment  
- **Proximal binding counts**  
- **phastCons averages**  
It then orders clusters by a chosen metric, applies user‐specified deletions/flips, builds a continuous cluster mapping, and updates both `architectureDetails.txt` and the original BED.

---

### Inputs

| Flag                | Parameter         | Type    | Required | Description                                                                                       |
|---------------------|-------------------|---------|----------|---------------------------------------------------------------------------------------------------|
| `-t`, `--tss_bed`   | `<tss_bed>`       | string  | no       | Path to TSS BED file (for gene‐distance & hypergeom).                                             |
| `-n`, `--nplb_bed`  | `<nplb_bed>`      | string  | no       | Path to `nplb_clustered.bed` (for binding counts & mapping).                                      |
| `-p`, `--phastcons_bw` | `<bw>`         | string  | no       | Path to phastCons BigWig (for conservation averages).                                             |
| `-m`, `--metric`    | `{avg_dist,avg_binding,avg_phastcons,hybrid,none}` | string  | no (default=none) | Which metric to use for ordering clusters. “hybrid” combines normalized distance & binding.  |
| `-d`, `--delete_clusters` | `<list>`   | string  | no       | Comma‐separated original cluster IDs to drop (mapped to `None`).                                   |
| `-f`, `--flip_clusters`   | `<list>`   | string  | no       | Comma‐separated original cluster IDs whose strand should be flipped in the updated BED.           |
| `--cluster_map_tsv` | `<map_tsv>`       | string  | no       | Optional TSV mapping original→mapped IDs with `flip`; requires `-m none` if provided.              |

---

### Pipeline Steps

1. **Gene‐TSS Distance & Hypergeom**  
   Runs if both `--tss_bed` and `--nplb_bed` are given:  
   - Sort TSS (`preprocess_tss`)  
   - Compute closest gene distances (`compute_closest`)  
   - Filter columns, run hypergeom tests, plot p-value heatmaps  
   - Compute `avg_dist` per cluster  

2. **Proximal Binding Count**  
   Runs if `--nplb_bed` is given:  
   - Extend windows by 20 kb (`extend_cluster_windows`)  
   - Count overlaps (`count_proximal_binding`)  
   - Compute `avg_binding` per cluster  

3. **phastCons Average**  
   Runs if `--phastcons_bw` is given:  
   - Compute mean phastCons via `pyBigWig` helper  
   - Compute `avg_phast` per cluster  

4. **Metric Mapping & Ordering**  
   - Build a DataFrame of clusters and their metrics  
   - Normalize and/or combine metrics if `hybrid`  
   - Sort clusters according to the chosen metric  

5. **Build & Save Cluster Mapping**  
   - Drop deleted clusters  
   - Apply strand flips  
   - Assign continuous `mapped_cluster_id` from 1…N  
   - Save `cluster_mapping.tsv`  

6. **Update Outputs**  
   - Use the provided or newly computed map to:  
     - Rewrite `architectureDetails_updated.txt` with integer IDs  
     - Rewrite `nplb_clustered_updated.bed` with mapped IDs and flipped strands  

---

### Usage Example

```bash
python nplb_ordering.py \
  --tss_bed /path/to/tss.bed \
  --nplb_bed /path/to/nplb_clustered.bed \
  --phastcons_bw /path/to/phastcons.bw \
  -m hybrid \
  -d 1,2 \
  -f 3,4 \
  --cluster_map_tsv /path/to/cluster_mapping.tsv
```

> *Notes:*  
> - Omit any of `-t`, `-n`, `-p` to skip that analysis step.  
> - If using `--cluster_map_tsv`, set `-m none`.

## 5. Conceptual & Usage Overview of `nplb_classify.py`

This script builds an NPLB clustered‐window BED file from `architectureDetails.txt`, then applies a user‐supplied cluster‐mapping TSV to:

1. **Update** `architectureDetails.txt` → `architectureDetails_updated.txt` (integer IDs).  
2. **Generate** `nplb_clustered.bed` from the raw architecture details.  
3. **Remap & flip** the BED into `nplb_clustered_classified.bed` using the mapping TSV.

---

### Inputs

| Flag                              | Parameter              | Type   | Required | Description                                                                                      |
|-----------------------------------|------------------------|--------|----------|--------------------------------------------------------------------------------------------------|
| `-o`, `--nplb_classify_dir`       | `<classify_dir>`       | string | yes      | Directory containing `architectureDetails.txt` (output of promoterLearn).                       |
| `--cluster_map_tsv`               | `<map_tsv>`            | string | no       | TSV with columns `original_cluster_id`, `mapped_cluster_id`, `flip`.                            |

---

### Processing Steps

1. **Build `nplb_clustered.bed`**  
   - Read `architectureDetails.txt` from `nplb_classify_dir`.  
   - Parse each line (`chr:start-end`) into `chr`, `start`, `end`, infer `strand` (+/–).  
   - Use column 1 as `cluster`.  
   - Write `nplb_clustered.bed`.

2. **Update Architecture Details**  
   - Call `update_architecture_details()` to map and cast cluster IDs to integers in `architectureDetails_updated.txt`.

3. **Remap & Flip BED**  
   - Call `update_nplb_bed()` to read `nplb_clustered.bed`, apply the TSV mapping (drop `None`), flip strands where indicated, and write `nplb_clustered_classified.bed`.

---

### Usage

```bash
python nplb_classify.py \
  --nplb_classify_dir /path/to/NPLB/output \
  --cluster_map_tsv /path/to/cluster_mapping.tsv
```
