# Loads and joins the two AMD data sources (RNA-seq CPM expression + clinical/
# genotype metadata) into the numpy arrays the xvae/DEC pipeline trains on.
# This is the only file that touches the raw CSVs - run_xdec.py calls
# get_data() from here rather than reading files directly.
import os
import re

import numpy as np
import pandas as pd
from sklearn.preprocessing import LabelEncoder

_HERE = os.path.dirname(os.path.abspath(__file__))
# aak100_cpmdat.csv's first column has no header ("Unnamed: 0") rather than a
# literal 'sample_id' column - _read_rna() below renames it before indexing.
DEFAULT_RNA_FILE = os.path.join(_HERE, '..', '..', 'Bayesian Hyperparameter Tuning',
                                'data', 'aak100_cpmdat.csv')
DEFAULT_CLIN_FILE = os.path.join(_HERE, '..', '..', 'Dataset', 'MetaSheet_1_4.csv')

# Columns considered "core" clinical/genotype signal - kept regardless of prevalence.
CORE_CLINICAL_COLS = [
    'age', 'sex_F', 'sex_M',
    'A69S_GG', 'A69S_GT', 'A69S_TT',
    'Y402H_CC', 'Y402H_CT', 'Y402H_TT',
]

# Free-text-derived oc_/mh_ history flags mentioning AMD/macular degeneration are a
# near-direct proxy for the mgs_level label and would leak the outcome into the
# clinical features, so they're excluded from the 'curated' clinical feature set.
DEFAULT_LEAKAGE_REGEX = r'amd|macular degener'


def _read_rna(rna_file):
    """aak100_cpmdat.csv ships with a blank/'Unnamed: 0' first column holding
    the sample id, not a literally-named 'sample_id' column. Normalise that
    before indexing so downstream code can always rely on 'sample_id'."""
    rna_df = pd.read_csv(rna_file)
    if 'sample_id' not in rna_df.columns:
        rna_df = rna_df.rename(columns={rna_df.columns[0]: 'sample_id'})
    return rna_df


def _build_curated_clinical(clin_df, min_history_prevalence=0.05,
                             leakage_regex=DEFAULT_LEAKAGE_REGEX):
    """Age + sex + genotype dummies, plus prevalence-filtered oc_/mh_ history flags
    with AMD-mentioning columns excluded as label leakage.
    """
    core_cols = [c for c in CORE_CLINICAL_COLS if c in clin_df.columns]
    history_cols = [c for c in clin_df.columns if c.startswith('oc_') or c.startswith('mh_')]

    prevalence = clin_df[history_cols].mean(axis=0)
    leak = re.compile(leakage_regex, re.IGNORECASE)
    kept_history = [
        c for c in history_cols
        if prevalence[c] >= min_history_prevalence and not leak.search(c)
    ]

    clin = clin_df[core_cols + kept_history].astype(np.float32)
    if 'age' in clin.columns:
        # Everything else here is a 0/1 dummy (sex, genotype, oc_/mh_ flags), but
        # age is continuous - min-max scale it into [0, 1] so it fits the binary/
        # sigmoid+BCE clinical branch's assumptions instead of corrupting it.
        age_range = clin['age'].max() - clin['age'].min()
        clin['age'] = (clin['age'] - clin['age'].min()) / age_range if age_range > 0 else 0.0
    return clin


def get_data(rna_file=DEFAULT_RNA_FILE, clin_file=DEFAULT_CLIN_FILE, label_col='mgs_level',
             clinical_mode='curated', min_history_prevalence=0.05,
             leakage_regex=DEFAULT_LEAKAGE_REGEX):
    """Load and join the AMD RNA-seq CPM matrix with the matched clinical sheet.

    Both files are indexed by `sample_id`; only samples present in both are kept.

    clinical_mode:
        'full'    - use every clinical column as-is (already one-hot/dummy encoded).
                    MetaSheet_1_4.csv has ~600+ mostly free-text-derived history
                    columns, most near-empty - for n~150 samples this makes the
                    clinical block far higher-dimensional than the sample count.
        'curated' - age/sex/genotype + prevalence-filtered (>= min_history_prevalence)
                    oc_/mh_ history flags, excluding any matching `leakage_regex`
                    (default excludes AMD/macular-degeneration mentions, which would
                    leak the label). Recommended default for clustering on this cohort.
    """
    rna_df = _read_rna(rna_file).set_index('sample_id')
    clin_df = pd.read_csv(clin_file).set_index('sample_id')

    # Only samples present in BOTH files can be used.
    common_ids = rna_df.index.intersection(clin_df.index)
    if len(common_ids) == 0:
        raise ValueError(
            "No overlapping sample_id between '{}' and '{}'".format(rna_file, clin_file)
        )
    rna_df = rna_df.loc[common_ids]
    clin_df = clin_df.loc[common_ids]

    # The label (mgs_level) lives in the RNA file; strip it from both frames
    # so it never ends up as a model input feature by accident. `y` is kept
    # only for reporting ARI/NMI/accuracy against DEC's unsupervised clusters
    # afterward - DEC never sees these labels during training.
    labels_raw = rna_df[label_col].astype(str)
    rna = rna_df.drop(columns=[label_col]).astype(np.float32)
    clin_drop = [c for c in [label_col] if c in clin_df.columns]
    clin_df = clin_df.drop(columns=clin_drop)

    if clinical_mode == 'curated':
        clin = _build_curated_clinical(clin_df, min_history_prevalence, leakage_regex)
    elif clinical_mode == 'full':
        clin = clin_df.astype(np.float32)
    else:
        raise ValueError("clinical_mode must be 'full' or 'curated', got {!r}".format(clinical_mode))

    le = LabelEncoder()
    y = le.fit_transform(labels_raw)

    return {
        'sample_id': common_ids.to_numpy(),
        'rnanp': rna.to_numpy(),
        'clin': clin.to_numpy(),
        'y': y,
        'label_classes': le.classes_,
        'gene_cols': list(rna.columns),
        'clin_cols': list(clin.columns),
    }
