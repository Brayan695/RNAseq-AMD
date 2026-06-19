"""
Linear Discriminant Analysis (LDA) for AMD Classification
==========================================================
Dataset : aak100_cpmdat.csv  (CPM-normalised RNA-seq, 81 genes, 166 samples)
Classes : MGS1 (early/normal, n=105) vs MGS4 (advanced AMD, n=61)

Method reference
----------------
Fisher, R.A. (1936). "The use of multiple measurements in taxonomic problems."
See also: https://en.wikipedia.org/wiki/Linear_discriminant_analysis

LDA finds the linear combination of features w that maximises the ratio
of between-class variance to within-class variance (Fisher's criterion):

    J(w) = (w^T S_B w) / (w^T S_W w)

Number of components
--------------------
The between-class scatter S_B is a sum of n_classes rank-<=1 matrices, so
rank(S_B) <= n_classes - 1.  This is a binary problem (MGS1 vs MGS4), hence
n_classes - 1 = 1: there is exactly ONE discriminant direction.  No package
(sklearn, R's MASS, mlxtend, ...) can produce a meaningful second LDA
component on 2-class data, so we use a single component throughout.
(sklearn user guide: n_components <= min(n_classes - 1, n_features).)

Implementation uses sklearn's LinearDiscriminantAnalysis (as in the linked
GeeksforGeeks tutorial); classification uses a NearestCentroid on the 1-D
projection; ROC / AUC use the signed Component-1 score.

Outputs
-------
lda_nc1.png          -- ROC | Comp1 dist | Confusion matrix | Gene coefficients
lda_coefficients.csv -- discriminant coefficients (Component 1)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.neighbors import NearestCentroid
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

COLORS = {"MGS1": "#1f77b4", "MGS4": "#ff7f0e"}


# ══════════════════════════════════════════════════════════════════════════════
# Single-component LDA  (sklearn LDA + NearestCentroid)
# ══════════════════════════════════════════════════════════════════════════════
class LDA1:
    """
    Single-component LDA built on sklearn's LinearDiscriminantAnalysis.

    Binary problem -> exactly one Fisher discriminant (n_classes - 1 = 1).
    Classification uses NearestCentroid on the 1-D projection.
    ROC / decision score uses the signed Component-1 projection.
    """

    def fit(self, X, y):
        self.lda_ = LinearDiscriminantAnalysis(n_components=1)
        z1 = self.lda_.fit_transform(X, y)[:, 0]            # (n,)

        # Sign-correct: ensure class 1 (MGS4) has the higher LD1 mean
        self.sign_ = 1.0 if z1[y == 1].mean() >= z1[y == 0].mean() else -1.0
        z1 = z1 * self.sign_

        self.clf_  = NearestCentroid().fit(z1[:, None], y)
        self.coef_ = self.lda_.scalings_[:, 0] * self.sign_  # LD1 loadings
        return self

    def decision_score(self, X):
        """Signed Component-1 score (used for ROC)."""
        return self.lda_.transform(X)[:, 0] * self.sign_

    def predict(self, X):
        return self.clf_.predict(self.decision_score(X)[:, None])


# ══════════════════════════════════════════════════════════════════════════════
# 1. Load & scale
# ══════════════════════════════════════════════════════════════════════════════
df = pd.read_csv(DATA_PATH, index_col=0)
feature_cols = [c for c in df.columns if c != "mgs_level"]
X = df[feature_cols].values
y = df["mgs_level"].values

le       = LabelEncoder()
y_enc    = le.fit_transform(y)            # MGS1->0, MGS4->1
scaler   = StandardScaler()
X_scaled = scaler.fit_transform(X)

print("=" * 60)
print("Linear Discriminant Analysis -- AMD (MGS1 vs MGS4)")
print("=" * 60)
print(f"Classes  : {le.classes_}")
print(f"Samples  : {len(y_enc)}  (MGS1={(y_enc==0).sum()}, MGS4={(y_enc==1).sum()})")
print(f"Features : {X.shape[1]}")
print("Components: 1  (binary problem -> n_classes - 1 = 1)")


# ══════════════════════════════════════════════════════════════════════════════
# 2. Cross-validation
# ══════════════════════════════════════════════════════════════════════════════
def run_cv_lda(X, y, cv):
    """
    10-fold CV with single-component LDA1.

    Returns
    -------
    y_score : (n,)  Component-1 scores across all folds (for ROC)
    y_pred  : (n,)  NearestCentroid predictions
    """
    n       = len(y)
    y_score = np.zeros(n)
    y_pred  = np.zeros(n, dtype=int)

    for tr, te in cv.split(X, y):
        clf = LDA1().fit(X[tr], y[tr])
        y_score[te] = clf.decision_score(X[te])
        y_pred[te]  = clf.predict(X[te])

    return y_score, y_pred


cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

y_score, y_pred = run_cv_lda(X_scaled, y_enc, cv)
fpr, tpr, _     = roc_curve(y_enc, y_score)
roc_auc         = auc(fpr, tpr)

print(f"\nROC-AUC = {roc_auc:.4f}  [10-fold CV]")
print(classification_report(y_enc, y_pred, target_names=le.classes_))


# ══════════════════════════════════════════════════════════════════════════════
# 3. Fit on full data for projection / coefficients
# ══════════════════════════════════════════════════════════════════════════════
lda_full = LDA1().fit(X_scaled, y_enc)
Z1       = lda_full.decision_score(X_scaled)        # (n,) Component-1 score

coef_df = pd.DataFrame({
    "Gene":        feature_cols,
    "Coefficient": lda_full.coef_,
}).sort_values("Coefficient", key=abs, ascending=False).reset_index(drop=True)

print("\nTop 20 discriminant genes (|Component 1 coefficient|):")
print(coef_df.head(20).to_string(index=False))

coef_df.to_csv(f"{OUT_DIR}\\lda_coefficients.csv", index=False)
print("\nCoefficients saved -> lda_coefficients.csv")


# ══════════════════════════════════════════════════════════════════════════════
# 4. Generate the 4-panel figure
# ══════════════════════════════════════════════════════════════════════════════
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
ax_roc = axes[0, 0];  ax_dist = axes[0, 1]
ax_cm  = axes[1, 0];  ax_d    = axes[1, 1]

fig.suptitle(
    "LDA -- AMD (MGS1 vs MGS4)  |  n_components = 1\n"
    f"ROC-AUC = {roc_auc:.3f}  [10-fold CV]",
    fontsize=13, fontweight="bold",
)

# ── Panel A: ROC ────────────────────────────────────────────────────────────
ax_roc.plot(fpr, tpr, color="steelblue", lw=2.5,
            label=f"LDA  (AUC = {roc_auc:.3f})")
ax_roc.fill_between(fpr, tpr, alpha=0.08, color="steelblue")
ax_roc.plot([0, 1], [0, 1], "k--", lw=1)
ax_roc.set_xlabel("False Positive Rate")
ax_roc.set_ylabel("True Positive Rate")
ax_roc.set_title("(A)  ROC Curve  [10-fold CV]")
ax_roc.legend(loc="lower right", fontsize=10)
ax_roc.grid(alpha=0.3)

# ── Panel B: Component 1 score distribution ─────────────────────────────────
for cls_idx, cls_name in enumerate(le.classes_):
    ax_dist.hist(Z1[y_enc == cls_idx], bins=20, alpha=0.55,
                 density=True, color=COLORS[cls_name], label=cls_name,
                 edgecolor="white", linewidth=0.4)
ax_dist.set_xlabel("Component 1 score")
ax_dist.set_ylabel("Density")
ax_dist.set_title("(B)  Component 1 Score Distributions")
ax_dist.legend();  ax_dist.grid(alpha=0.3)

# ── Panel C: Confusion matrix ───────────────────────────────────────────────
cm   = confusion_matrix(y_enc, y_pred)
disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=le.classes_)
disp.plot(ax=ax_cm, colorbar=False, cmap="Blues")
ax_cm.set_title("(C)  Confusion Matrix  [10-fold CV]")

# ── Panel D: Gene coefficients ──────────────────────────────────────────────
top20  = coef_df.head(20).copy()
colors = ["#d62728" if v > 0 else "#1f77b4" for v in top20["Coefficient"]]
ax_d.barh(top20["Gene"][::-1], top20["Coefficient"][::-1],
          color=colors[::-1], edgecolor="white", linewidth=0.4)
ax_d.axvline(0, color="black", linewidth=0.8)
ax_d.set_xlabel("Discriminant Coefficient (standardised)")
ax_d.set_title("(D)  Top 20 Genes by |Coefficient|\n"
               "(red = higher in MGS4,  blue = higher in MGS1)")
ax_d.grid(axis="x", alpha=0.3)

plt.tight_layout()
out_path = f"{OUT_DIR}\\lda_nc1.png"
fig.savefig(out_path, dpi=150, bbox_inches="tight")
plt.show()
print(f"\nFigure saved -> lda_nc1.png")
print("\nDone.")
