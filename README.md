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
## 1. Conceptual & Usage Overview of `nplb_train.py`

This Python script parses the NPLB `promoterLearn` training output (`architectureDetails.txt`) and converts it into a BED file of clustered windows.

### Inputs

| Flag                           | Parameter           | Type   | Required | Description                                                                                               |
|--------------------------------|---------------------|--------|----------|-----------------------------------------------------------------------------------------------------------|
| `--output_prefix`, `-o`        | `<output_prefix>`   | string | yes      | Output prefix (directory + basename) where `architectureDetails.txt` was written by `promoterLearn`.     |

### Functions Utilized

- `parse_nplb_train_args()`: parses the `--output_prefix` argument.  
  *(defined in `nplb_helper_functions.py`)*

### Processing Steps

1. **Locate** `architectureDetails.txt` in the directory of `output_prefix`.  
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
  --output_prefix /home/samarth/ctcf_sequence_data/output/NPLB/HFF/output_prefix
```
