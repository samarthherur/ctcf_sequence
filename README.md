## Conceptual Overview of `fimo_neighbourhood_analysis.py`

This script implements an end-to-end pipeline to take raw FIMO motif calls and integrate them with ChIP-seq and boundary data, producing statistical summaries and visualizations at each step.

### 1. Inputs
- **FIMO output** (TSV of motif hits)
- **Genome FASTA** (reference sequence)
- **ChIP-seq peaks** (BED file)
- **Boundary regions** (BED file)
- **Neighborhood size** (window around each motif, e.g. 50 bp)
- **Output directory** (where all results are saved)

---

### 2. Neighborhood Extraction
1. **Extend** each FIMO hit by ±_neighborhood_size_ around its genomic coordinates.
2. **Merge** overlapping windows so each motif has a unique “neighborhood” window.
3. **Split** by strand and write plus/minus BED files.
4. **Fetch** the corresponding sequences with `bedtools getfasta`.
5. **Save** a combined FASTA of all neighborhoods.

---

### 3. Visualization & Initial Statistics
- **Heatmap & Logo** of the raw neighborhood FASTA.
- **Median score** of all neighborhood windows.

---

### 4. ChIP-seq Filtering
- **Intersect** neighborhood windows with ChIP-seq peaks.
- **Save** filtered BED & FASTA.
- **Redo** heatmap & logo on filtered sequences.

---

### 5. Boundary Intersection
- **Split** filtered windows into:
  - Those that **intersect** boundary regions  
  - Those that **do not** intersect
- **Heatmap & logo** for both subsets.

---

### 6. Score Distributions & Tests
- **Compute** median scores for intersecting vs non-intersecting windows.
- **Plot** and save:
  - Density overlay (`score_distribution.png`)
  - Overlaid histograms (`score_histogram.png`)
- **Kolmogorov–Smirnov test** between the two groups.

---

### 7. Motif Classification
- **Filter** the original FIMO TSV by ChIP peaks.
- **Assign** each motif to boundary windows.
- **Classify** motifs into:
  1. Class 1 (outside any boundary)  
  2. Class 2 (inside a boundary, no special neighbor relationship)  
  3. Class 3 (strand-paired motifs across adjacent boundaries)
- **Write** `classified_motifs.tsv`.

---

### 8. Classified Motif Summaries
- **Total count** and **per-class counts**.
- **Median score** per class.
- **KS tests** between every pair of classes.
- **Density plot** of class-specific score distributions.

---

### 9. Repeat-Filtering & Final Visuals
- **Remove** any sequences with `N`/`n` or lowercase letters.
- **Report** number of repeat-free sequences and per-base counts.
- **Save** a `_wr.fasta` of cleaned sequences.
- **Heatmap & logo** of repeat-free FASTA.

---

### 10. Neighborhood Matrix
- **Count** overlaps with ChIP and boundary for each window (`-c`).
- **Save** `neighbourhood_matrix.tsv` for downstream analyses.

---

_All results (metrics, statistics, and plots) are logged to `summary.txt` and accompanying PNG files in the specified output folder._  