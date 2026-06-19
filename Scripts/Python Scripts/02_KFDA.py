"""
Kernel Fisher Discriminant Analysis (KFDA) for AMD Classification
==================================================================
Dataset : aak100_cpmdat.csv  (CPM-normalised RNA-seq, 81 genes, 166 samples)
Classes : MGS1 (early/normal, n=105) vs MGS4 (advanced AMD, n=61)

Method reference
----------------
Mika et al. (1999). "Fisher Discriminant Analysis with Kernels."
See also: https://en.wikipedia.org/wiki/Kernel_Fisher_discriminant_analysis

Package used
------------
    pip install kfda

Number of components
--------------------
KFDA is Fisher discriminant analysis performed in a kernel feature space, so
it obeys the SAME rank limit as ordinary LDA: the between-class scatter has
rank <= n_classes - 1.  This is a binary problem (MGS1 vs MGS4), hence
n_classes - 1 = 1 -- exactly ONE meaningful discriminant direction.  The kfda
package will numerically return extra eigenvectors if asked, but for 2-class
data they are near-zero-eigenvalue noise, not discriminants.  We therefore use
a single component and compare all four kernels at n_components = 1.

ROC / AUC  : driven by Component 1 (the discriminant direction)
Predictions: NearestCentroid on Component 1 (built into kfda)

Outputs
-------
kfda_nc1.png  -- ROC all kernels | Comp1 dist | Confusion matrix | AUC bar
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

from kfda import Kfda
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import (
    roc_curve, auc, classification_report,
    confusion_matrix, ConfusionMatrixDisplay,
)

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_PATH = (
    r"C:\Users\Brayan Gutierrez\Desktop\RNAseq-AMD\Dataset\aak100_cpmdat.csv"
)
OUT_DIR = r"C:\Users\Brayan Gutierrez\Desktop\RNAseq-AMD\Scripts\Python Scripts"

COLORS  = {"MGS1": "#1f77b4", "MGS4": "#ff7f0e"}
PALETTE = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]

# ── Kernel configurations (single component) ──────────────────────────────────
N_COMPONENTS = 1
KERNEL_DEFS = {
    "RBF (g=0.01)":     {"kernel": "rbf",    "gamma": 0.01},
    "RBF (g=0.001)":    {"kernel": "rbf",    "gamma": 0.001},
    "Polynomial (d=2)": {"kernel": "poly",   "degree": 2},
    "Linear":           {"kernel": "linear"},
}

# ══════════════════════════════════════════════════════════════════════════════
# 1. Load & scale
# ══════════════════════════════════════════════════════════════════════════════
df = pd.read_csv(DATA_PATH, index_col=0)
feature_cols = [c for c in df.columns if c != "mgs_level"]
X = df[feature_cols].values
y = df["mgs_level"].values

le       = LabelEncoder()
y_enc    = le.fit_transform(y)
scaler   = StandardScaler()
X_scaled = scaler.fit_transform(X)

print("=" * 60)
print("Kernel Fisher Discriminant Analysis -- AMD (MGS1 vs MGS4)")
print("=" * 60)
print(f"Classes  : {le.classes_}")
print(f"Samples  : {len(y_enc)}  (MGS1={(y_enc==0).sum()}, MGS4={(y_enc==1).sum()})")
print(f"Features : {X.shape[1]}")
print("Components: 1  (binary problem -> n_classes - 1 = 1)")


# ══════════════════════════════════════════════════════════════════════════════
# 2. Cross-validation helper
# ══════════════════════════════════════════════════════════════════════════════
def run_cv(kernel_kwargs, X, y, cv):
    """
    10-fold CV for one Kfda configuration.

    Returns
    -------
    y_score : Component-1 projection (sign-corrected; used for ROC)
    y_pred  : NearestCentroid predictions
    """
    n       = len(y)
    y_score = np.zeros(n)
    y_pred  = np.zeros(n, dtype=int)

    for tr, te in cv.split(X, y):
        clf  = Kfda(**kernel_kwargs)
        clf.fit(X[tr], y[tr])

        Z_tr = clf.transform(X[tr])
        Z_te = clf.transform(X[te])

        # Sign-correct component 1 so MGS4 > MGS1
        c0 = Z_tr[y[tr] == 0, 0].mean()
        c1 = Z_tr[y[tr] == 1, 0].mean()
        sign = np.sign(c1 - c0) if c1 != c0 else 1.0

        y_score[te] = sign * Z_te[:, 0]
        y_pred[te]  = clf.predict(X[te])

    return y_score, y_pred


def fit_full(kernel_kwargs, X, y):
    """Fit on full data; sign-correct Component 1 for plotting."""
    clf = Kfda(**kernel_kwargs)
    clf.fit(X, y)
    z = clf.transform(X)[:, 0]
    if z[y == 1].mean() < z[y == 0].mean():
        z = -z
    return z


cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)


# ══════════════════════════════════════════════════════════════════════════════
# 3. Run all kernels at n_components = 1
# ══════════════════════════════════════════════════════════════════════════════
print(f"\n{'='*60}")
print(f"n_components = {N_COMPONENTS}")
print(f"{'='*60}")

kernel_results = {}
for name, kdef in KERNEL_DEFS.items():
    kw              = {**kdef, "n_components": N_COMPONENTS}
    y_score, y_pred = run_cv(kw, X_scaled, y_enc, cv)
    fpr, tpr, _     = roc_curve(y_enc, y_score)
    roc_auc         = auc(fpr, tpr)
    kernel_results[name] = dict(fpr=fpr, tpr=tpr, auc=roc_auc,
                                y_score=y_score, y_pred=y_pred, kw=kw)
    print(f"  {name:<22s}  AUC = {roc_auc:.4f}")

best_name = max(kernel_results, key=lambda k: kernel_results[k]["auc"])
best      = kernel_results[best_name]

print(f"\nBest kernel: {best_name}  (AUC = {best['auc']:.4f})")
print(f"\nClassification Report  [{best_name}]:")
print(classification_report(y_enc, best["y_pred"], target_names=le.classes_))

# Full-data Component-1 projection with best kernel (for distribution panel)
Z1 = fit_full(best["kw"], X_scaled, y_enc)


# ══════════════════════════════════════════════════════════════════════════════
# 4. Generate the 4-panel figure
# ══════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
ax_roc = axes[0, 0];  ax_dist = axes[0, 1]
ax_cm  = axes[1, 0];  ax_d    = axes[1, 1]

fig.suptitle(
    "KFDA -- AMD (MGS1 vs MGS4)  |  n_components = 1\n"
    f"Best kernel: {best_name}  |  AUC = {best['auc']:.3f}  [10-fold CV]",
    fontsize=12, fontweight="bold",
)

# ── Panel A: ROC all kernels ────────────────────────────────────────────────
for (name, res), col in zip(kernel_results.items(), PALETTE):
    ax_roc.plot(res["fpr"], res["tpr"], lw=2, color=col,
                label=f"{name}  ({res['auc']:.3f})")
ax_roc.plot([0, 1], [0, 1], "k--", lw=1)
ax_roc.set_xlabel("False Positive Rate")
ax_roc.set_ylabel("True Positive Rate")
ax_roc.set_title("(A)  ROC Curves  [10-fold CV]")
ax_roc.legend(loc="lower right", fontsize=8)
ax_roc.grid(alpha=0.3)

# ── Panel B: Component 1 score distribution (best kernel) ───────────────────
for cls_idx, cls_name in enumerate(le.classes_):
    ax_dist.hist(Z1[y_enc == cls_idx], bins=20, alpha=0.55,
                 density=True, color=COLORS[cls_name], label=cls_name,
                 edgecolor="white", linewidth=0.4)
ax_dist.set_xlabel("Component 1 score")
ax_dist.set_ylabel("Density")
ax_dist.set_title(f"(B)  Component 1 Distributions\n[Best kernel: {best_name}]")
ax_dist.legend();  ax_dist.grid(alpha=0.3)

# ── Panel C: Confusion matrix (best kernel) ─────────────────────────────────
cm   = confusion_matrix(y_enc, best["y_pred"])
disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=le.classes_)
disp.plot(ax=ax_cm, colorbar=False, cmap="Blues")
ax_cm.set_title(f"(C)  Confusion Matrix  [10-fold CV, {best_name}]")

# ── Panel D: AUC comparison across kernels ──────────────────────────────────
names_list = list(kernel_results.keys())
aucs_list  = [kernel_results[n]["auc"] for n in names_list]
bars = ax_d.bar(range(len(names_list)), aucs_list,
                color=PALETTE, alpha=0.82, edgecolor="white")
ax_d.set_xticks(range(len(names_list)))
ax_d.set_xticklabels(names_list, rotation=22, ha="right", fontsize=9)
ax_d.set_ylabel("AUC");  ax_d.set_ylim(0.45, 1.05)
ax_d.set_title("(D)  AUC Comparison Across Kernels")
ax_d.axhline(0.5, color="gray", linestyle="--", linewidth=1)
for bar, val in zip(bars, aucs_list):
    ax_d.text(bar.get_x() + bar.get_width() / 2,
              bar.get_height() + 0.01,
              f"{val:.3f}", ha="center", va="bottom", fontsize=9)
ax_d.grid(axis="y", alpha=0.3)

plt.tight_layout()
out_path = f"{OUT_DIR}\\kfda_nc1.png"
fig.savefig(out_path, dpi=150, bbox_inches="tight")
plt.show()
print(f"\nFigure saved -> kfda_nc1.png")
print("\nDone.")
