# Loads data for the ALL-GENES 4-class variant of CNC-VAE: unlike CNC-VAE-4class
# (which pulls out just the same 81 genes as the original binary CNC-VAE),
# this uses every gene in gene_input.csv directly (~18,056 genes, all 453
# donors) - no gene-matching step needed since there's no smaller reference
# gene set to align to. This is a high-dimensional, low-sample regime, so
# treat the PCA/raw baselines in compare_representations.py as load-bearing.
import os
import re

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import LabelEncoder

_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_GENE_INPUT_FILE = os.path.join(_HERE, '..', '..', 'Dataset', 'gene_input.csv')
DEFAULT_CLIN_FILE = os.path.join(_HERE, '..', '..', 'Dataset', 'MetaSheet_Processed.csv')

# genotype/sex columns are already one-hot encoded in the clinical sheet.
CORE_CLINICAL_PREFIXES = ('A69S_', 'Y402H_', 'sex_')

# Free-text-derived oc_/mh_ history flags mentioning AMD/macular degeneration are a
# near-direct proxy for the mgs_level label and would leak the outcome into the
# clinical arm, so they're excluded from the 'curated' clinical feature set.
DEFAULT_LEAKAGE_REGEX = r'amd|macular degener'


def _build_curated_clinical(clin_df, min_history_prevalence=0.05, leakage_regex=DEFAULT_LEAKAGE_REGEX):
    """Same curated-clinical-feature logic as the other CNC-VAE variants: keep
    age/sex/genotype always, plus oc_/mh_ history flags common enough to be
    signal rather than noise, minus anything that would leak the AMD diagnosis."""
    core_cols = [c for c in clin_df.columns if c == 'age' or c.startswith(CORE_CLINICAL_PREFIXES)]
    history_cols = [c for c in clin_df.columns if c.startswith('oc_') or c.startswith('mh_')]

    prevalence = clin_df[history_cols].mean(axis=0)
    leak = re.compile(leakage_regex, re.IGNORECASE)
    kept_history = [
        c for c in history_cols
        if prevalence[c] >= min_history_prevalence and not leak.search(c)
    ]
    return clin_df[core_cols + kept_history].astype(np.float32)


def get_data(gene_input_file=DEFAULT_GENE_INPUT_FILE, clin_file=DEFAULT_CLIN_FILE,
             label_col='mgs_level', clinical_mode='curated', min_history_prevalence=0.05,
             leakage_regex=DEFAULT_LEAKAGE_REGEX):
    """Load the FULL gene_input.csv transcriptome (~18,056 genes, all 453 AMD
    donors), joined with clinical/genotype covariates, for the 4-class
    MGS1..MGS4 grade task.

    This is a high-dimensional, low-sample regime (~18k genes, 453 samples) -
    overfitting risk is high, so treat the PCA/raw baselines in
    compare_representations.py as load-bearing, not just a formality.
    """
    gi_df = pd.read_csv(gene_input_file).set_index('sample_id')
    clin_df = pd.read_csv(clin_file).set_index('sample_id')

    # Every column except the label is treated as a gene - no gene-matching
    # step needed here (contrast with CNC-VAE-4class, which has to align its
    # 81-gene subset against a differently-named reference file).
    gene_cols = [c for c in gi_df.columns if c != label_col]
    rna_df = gi_df[gene_cols]
    labels_df = gi_df[label_col]

    common_ids = rna_df.index.intersection(clin_df.index)
    if len(common_ids) == 0:
        raise ValueError(
            "No overlapping sample_id between '{}' and '{}'".format(gene_input_file, clin_file)
        )
    rna_df = rna_df.loc[common_ids].astype(np.float32)
    labels_raw = labels_df.loc[common_ids].astype(str)
    clin_df = clin_df.loc[common_ids]
    clin_drop = [c for c in [label_col] if c in clin_df.columns]
    clin_df = clin_df.drop(columns=clin_drop)

    if clinical_mode == 'curated':
        clin = _build_curated_clinical(clin_df, min_history_prevalence, leakage_regex)
    elif clinical_mode == 'full':
        clin = clin_df.astype(np.float32)
    else:
        raise ValueError("clinical_mode must be 'full' or 'curated', got {!r}".format(clinical_mode))

    # 4-class label: LabelEncoder sorts alphabetically, which happens to also be
    # severity order here (MGS1 < MGS2 < MGS3 < MGS4), so y=0..3 matches grade order.
    le = LabelEncoder()
    y = le.fit_transform(labels_raw)

    return {
        'sample_id': common_ids.to_numpy(),
        'rnanp': rna_df.to_numpy(),  # (453, ~18056) - much wider than the 81-gene variants
        'clin': clin.to_numpy(),
        'y': y,
        'label_classes': le.classes_,
        'gene_cols': list(rna_df.columns),
        'clin_cols': list(clin.columns),
    }


class AMDDataset:
    """Whole-cohort access plus stratified K-fold splits (stratified on `label_col`)
    for CNC-VAE Clin+mRNA integration on the full 453-donor, 4-class, all-genes
    AMD cohort.
    """

    def __init__(self, gene_input_file=DEFAULT_GENE_INPUT_FILE, clin_file=DEFAULT_CLIN_FILE,
                 label_col='mgs_level', n_splits=5, seed=42, clinical_mode='curated'):
        self.data = get_data(gene_input_file, clin_file, label_col, clinical_mode=clinical_mode)
        self.n_splits = n_splits
        self.seed = seed
        self._folds = None

    @property
    def whole(self):
        return self.data

    def fold(self, k):
        """1-indexed fold number. Returns (train_dict, test_dict)."""
        # Folds are computed once (on first call) and cached, so every call to
        # fold(k) during a run sees the exact same split.
        if self._folds is None:
            skf = StratifiedKFold(n_splits=self.n_splits, shuffle=True, random_state=self.seed)
            self._folds = list(skf.split(self.data['rnanp'], self.data['y']))
        train_idx, test_idx = self._folds[k - 1]
        return self._subset(train_idx), self._subset(test_idx)

    def _subset(self, idx):
        """Slice every array in self.data down to the given row indices."""
        d = self.data
        return {
            'sample_id': d['sample_id'][idx],
            'rnanp': d['rnanp'][idx],
            'clin': d['clin'][idx],
            'y': d['y'][idx],
        }
