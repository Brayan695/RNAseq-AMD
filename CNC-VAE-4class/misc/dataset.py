import os
import re

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import LabelEncoder

_HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_GENE_INPUT_FILE = os.path.join(_HERE, '..', '..', 'Dataset', 'gene_input.csv')
DEFAULT_AAK_FILE = os.path.join(_HERE, '..', '..', 'Dataset', 'aak100_cpmdat.csv')
DEFAULT_CLIN_FILE = os.path.join(_HERE, '..', '..', 'Dataset', 'MetaSheet_Processed.csv')

# genotype/sex columns are already one-hot encoded in the clinical sheet.
CORE_CLINICAL_PREFIXES = ('A69S_', 'Y402H_', 'sex_')

# Free-text-derived oc_/mh_ history flags mentioning AMD/macular degeneration are a
# near-direct proxy for the mgs_level label and would leak the outcome into the
# clinical arm, so they're excluded from the 'curated' clinical feature set.
DEFAULT_LEAKAGE_REGEX = r'amd|macular degener'


def _match_ensg_to_symbol(aak_df, gi_df, ens_order):
    """Identify which gene_input.csv SYMBOL column corresponds to each of the 81
    ENSG-named genes in aak100_cpmdat.csv, by matching expression VALUES on the
    166 samples the two files share (gene_input.csv has no direct ENSG<->symbol
    map, only gene symbols, and aak100_cpmdat.csv only covers the MGS1/MGS4 subset).
    """
    sym_cands = [c for c in gi_df.columns if c not in ('sample_id', 'mgs_level')]
    shared = aak_df[['sample_id'] + ens_order].merge(
        gi_df[['sample_id'] + sym_cands], on='sample_id'
    )
    if len(shared) == 0:
        raise ValueError("aak100_cpmdat.csv and gene_input.csv share no sample_id - can't match genes")

    A = shared[ens_order].to_numpy(float)
    G = shared[sym_cands].to_numpy(float)
    ens_to_sym = {}
    for j, e in enumerate(ens_order):
        d = np.abs(G - A[:, [j]]).max(axis=0)
        k = int(d.argmin())
        if d[k] > 1e-3:
            raise ValueError("No gene_input.csv column matches {} (best abs diff {:.4g})".format(e, d[k]))
        ens_to_sym[e] = sym_cands[k]
    if len(set(ens_to_sym.values())) != len(ens_order):
        raise ValueError("Two ENSG genes mapped to the same gene_input.csv symbol")
    return ens_to_sym


def _build_curated_clinical(clin_df, min_history_prevalence=0.05, leakage_regex=DEFAULT_LEAKAGE_REGEX):
    core_cols = [c for c in clin_df.columns if c == 'age' or c.startswith(CORE_CLINICAL_PREFIXES)]
    history_cols = [c for c in clin_df.columns if c.startswith('oc_') or c.startswith('mh_')]

    prevalence = clin_df[history_cols].mean(axis=0)
    leak = re.compile(leakage_regex, re.IGNORECASE)
    kept_history = [
        c for c in history_cols
        if prevalence[c] >= min_history_prevalence and not leak.search(c)
    ]
    return clin_df[core_cols + kept_history].astype(np.float32)


def get_data(gene_input_file=DEFAULT_GENE_INPUT_FILE, aak_file=DEFAULT_AAK_FILE,
             clin_file=DEFAULT_CLIN_FILE, label_col='mgs_level',
             clinical_mode='curated', min_history_prevalence=0.05,
             leakage_regex=DEFAULT_LEAKAGE_REGEX):
    """Load the 81-gene subset (matched from gene_input.csv against the ENSG ids in
    aak100_cpmdat.csv) for all 453 AMD donors, joined with clinical/genotype
    covariates, for the 4-class MGS1..MGS4 grade task.
    """
    gi_df = pd.read_csv(gene_input_file)
    aak_df = pd.read_csv(aak_file)
    clin_df = pd.read_csv(clin_file).set_index('sample_id')

    ens_order = [c for c in aak_df.columns if c.startswith('ENSG')]
    ens_to_sym = _match_ensg_to_symbol(aak_df, gi_df, ens_order)
    sym_order = [ens_to_sym[e] for e in ens_order]

    rna_df = gi_df[['sample_id'] + sym_order].copy()
    rna_df.columns = ['sample_id'] + ens_order
    rna_df = rna_df.set_index('sample_id')
    labels_df = gi_df.set_index('sample_id')[label_col]

    common_ids = rna_df.index.intersection(clin_df.index)
    if len(common_ids) == 0:
        raise ValueError(
            "No overlapping sample_id between gene_input.csv and '{}'".format(clin_file)
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

    le = LabelEncoder()
    y = le.fit_transform(labels_raw)

    return {
        'sample_id': common_ids.to_numpy(),
        'rnanp': rna_df.to_numpy(),
        'clin': clin.to_numpy(),
        'y': y,
        'label_classes': le.classes_,
        'gene_cols': list(rna_df.columns),
        'clin_cols': list(clin.columns),
    }


class AMDDataset:
    """Whole-cohort access plus stratified K-fold splits (stratified on `label_col`)
    for CNC-VAE Clin+mRNA integration on the full 453-donor, 4-class AMD cohort.
    """

    def __init__(self, gene_input_file=DEFAULT_GENE_INPUT_FILE, aak_file=DEFAULT_AAK_FILE,
                 clin_file=DEFAULT_CLIN_FILE, label_col='mgs_level', n_splits=5, seed=42,
                 clinical_mode='curated'):
        self.data = get_data(gene_input_file, aak_file, clin_file, label_col,
                              clinical_mode=clinical_mode)
        self.n_splits = n_splits
        self.seed = seed
        self._folds = None

    @property
    def whole(self):
        return self.data

    def fold(self, k):
        """1-indexed fold number. Returns (train_dict, test_dict)."""
        if self._folds is None:
            skf = StratifiedKFold(n_splits=self.n_splits, shuffle=True, random_state=self.seed)
            self._folds = list(skf.split(self.data['rnanp'], self.data['y']))
        train_idx, test_idx = self._folds[k - 1]
        return self._subset(train_idx), self._subset(test_idx)

    def _subset(self, idx):
        d = self.data
        return {
            'sample_id': d['sample_id'][idx],
            'rnanp': d['rnanp'][idx],
            'clin': d['clin'][idx],
            'y': d['y'][idx],
        }
