# Loads and joins the two AMD data sources (81-gene RNA-seq CPM expression +
# clinical/genotype metadata) into the numpy arrays the CNC-VAE model trains on.
# This is the ONLY file that touches the raw CSVs - every other script in this
# project (run_cncvae.py, compare_representations.py, analyse_representations.py)
# calls get_data()/AMDDataset from here rather than reading files directly.
import os
import re

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import LabelEncoder

# File paths resolve relative to this file's location, two directories up into a
# sibling "Dataset" folder (i.e. ../../Dataset relative to misc/), so they work
# regardless of what directory the calling script is run from.
_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_RNA_FILE = os.path.join(_HERE, '..', '..', '..', 'Dataset', 'aak100_cpmdat.csv')
DEFAULT_CLIN_FILE = os.path.join(_HERE, '..', '..', '..', 'Dataset', 'MetaSheet_1_4.csv')

# Prefixes for "core" clinical/genotype signal - kept regardless of prevalence.
# Prefix-matched (not an exact-name list) so it survives renames like sex_F/sex_M
# being consolidated into a single sex_bin column.
CORE_CLINICAL_PREFIXES = ('sex_', 'A69S_', 'Y402H_')

# Free-text-derived oc_/mh_ history flags mentioning AMD/macular degeneration are a
# near-direct proxy for the mgs_level label and would leak the outcome into the
# clinical arm, so they're excluded from the 'curated' clinical feature set.
DEFAULT_LEAKAGE_REGEX = r'amd|macular degener'


def _build_curated_clinical(clin_df, min_history_prevalence=0.05,
                             leakage_regex=DEFAULT_LEAKAGE_REGEX):
    """Age + sex + genotype dummies, plus prevalence-filtered oc_/mh_ history flags
    with AMD-mentioning columns excluded as label leakage.
    """
    # "Core" columns (age, sex, genotype) are always kept regardless of how common
    # they are - use `c == 'age' or c.startswith(...)` rather than an exact-name
    # list, so this survives future column renames (e.g. sex_F/sex_M -> sex_bin).
    core_cols = [c for c in clin_df.columns if c == 'age' or c.startswith(CORE_CLINICAL_PREFIXES)]

    # The clinical sheet also has hundreds of free-text-derived oc_/mh_ columns
    # (ocular/medical history), most of which are 0 for nearly every patient -
    # e.g. a flag that's only true for 1 out of 166 samples is just noise to a
    # model this small, so only keep the ones common enough to carry signal.
    history_cols = [c for c in clin_df.columns if c.startswith('oc_') or c.startswith('mh_')]
    prevalence = clin_df[history_cols].mean(axis=0)  # fraction of samples where the flag is 1

    # Some of those history flags literally say "AMD" or "macular degeneration" -
    # keeping those would let the model see a near-restatement of the label it's
    # supposed to predict, so they're dropped regardless of prevalence.
    leak = re.compile(leakage_regex, re.IGNORECASE)
    kept_history = [
        c for c in history_cols
        if prevalence[c] >= min_history_prevalence and not leak.search(c)
    ]

    clin = clin_df[core_cols + kept_history].astype(np.float32)
    return clin


def get_data(rna_file=DEFAULT_RNA_FILE, clin_file=DEFAULT_CLIN_FILE, label_col='mgs_level',
             clinical_mode='full', min_history_prevalence=0.05,
             leakage_regex=DEFAULT_LEAKAGE_REGEX):
    """Load and join the AMD RNA-seq CPM matrix with the matched clinical sheet.

    Both files are indexed by `sample_id`; only samples present in both are kept.

    clinical_mode:
        'full'    - use every clinical column as-is (already one-hot/dummy encoded).
        'curated' - age/sex/genotype + prevalence-filtered (>= min_history_prevalence)
                    oc_/mh_ history flags, excluding any matching `leakage_regex`
                    (default excludes AMD/macular-degeneration mentions, which would
                    leak the label). Mirrors the reference CNC_VAE_AMD_V3 notebook.
    """
    rna_df = pd.read_csv(rna_file).set_index('sample_id')
    clin_df = pd.read_csv(clin_file).set_index('sample_id')

    # Only samples present in BOTH files can be used (the expression file may
    # cover fewer donors than the clinical sheet, e.g. aak100_cpmdat.csv only has
    # the MGS1/MGS4 subset while MetaSheet_Processed.csv covers all 4 grades).
    common_ids = rna_df.index.intersection(clin_df.index)
    if len(common_ids) == 0:
        raise ValueError(
            "No overlapping sample_id between '{}' and '{}'".format(rna_file, clin_file)
        )
    rna_df = rna_df.loc[common_ids]
    clin_df = clin_df.loc[common_ids]

    # The label (mgs_level: MGS1/MGS4 etc.) lives in the RNA file; strip it out of
    # both frames so it never accidentally ends up as a model input feature.
    labels_raw = rna_df[label_col].astype(str)
    rna = rna_df.drop(columns=[label_col]).astype(np.float32)
    clin_drop = [c for c in [label_col] if c in clin_df.columns]
    clin_df = clin_df.drop(columns=clin_drop)

    if clinical_mode == 'curated':
        clin = _build_curated_clinical(clin_df, min_history_prevalence, leakage_regex)
    elif clinical_mode == 'full':
        clin = clin_df.astype(np.float32)  # every remaining clinical column, unfiltered
    else:
        raise ValueError("clinical_mode must be 'full' or 'curated', got {!r}".format(clinical_mode))

    # Turn string labels ('MGS1', 'MGS4', ...) into integers (0, 1, ...) for the
    # model; label_classes remembers the mapping so it can be reversed later
    # (e.g. to print human-readable class names in plots/reports).
    le = LabelEncoder()
    y = le.fit_transform(labels_raw)

    return {
        'sample_id': common_ids.to_numpy(),
        'rnanp': rna.to_numpy(),      # (n_samples, n_genes) gene expression matrix
        'clin': clin.to_numpy(),      # (n_samples, n_clinical_features) clinical matrix
        'y': y,                       # (n_samples,) integer-encoded labels
        'label_classes': le.classes_,  # y=0 means label_classes[0], etc.
        'gene_cols': list(rna.columns),
        'clin_cols': list(clin.columns),
    }


class AMDDataset:
    """Whole-dataset access plus stratified K-fold splits (stratified on `label_col`)
    for CNC-VAE Clin+mRNA integration on the AMD cohort.

    The original IntegrativeVAEs repo shipped pre-built METABRIC fold CSVs; we
    don't have those for this cohort, so folds are built on the fly with
    StratifiedKFold to preserve the MGS1/MGS4 class balance in each split.
    """

    def __init__(self, rna_file=DEFAULT_RNA_FILE, clin_file=DEFAULT_CLIN_FILE,
                 label_col='mgs_level', n_splits=5, seed=42, clinical_mode='full'):
        self.data = get_data(rna_file, clin_file, label_col, clinical_mode=clinical_mode)
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
        train_idx, test_idx = self._folds[k - 1]  # k is 1-indexed; list is 0-indexed
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
