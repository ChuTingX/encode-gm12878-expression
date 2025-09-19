# Methods Summary (GM12878)

## Data
- **Genome build**: hg19
- **Gene models (GTF)**: GENCODE v19
- **CAGE (expression proxy)**: RIKEN CAGE bigWig for GM12878 (nuclear, plus-strand signal)
- **Predictors**: ENCODE integration bigWigs (histone/DNase) for GM12878

## Preprocessing
- Import transcripts; keep those with length ≥ 4100 bp on selected chromosomes.
- Derive 1-bp TSS per transcript; compute CAGE RPM in ±50bp window; choose the **max-CAGE TSS per gene**.
- Around each chosen TSS, build `nbins_up + nbins_down + 1` bins (here 40+40+1) with strand awareness.

## Feature selection (best-bin)
- Split genes: D1 for selecting the **best bin per track**; D2 for model evaluation.
- For each predictor track, pick the bin maximizing |correlation| with expression using `log2(signal + 0.1)`.

## Modeling
- **Classification**: logistic regression, random forest, SVM → predict ON/OFF.
- **Regression**: LASSO, random forest, MARS, SVM on ON genes.
- **Two-step combos**: classifier gate → regressor on predicted-ON.
- **Evaluation**: 10-fold CV; report AUC/Misclassification for classification and RMSE/Pearson r for regression.

## Notes
- Large data are not versioned; see code under `R/` and `scripts/` and configs in `conf/`.