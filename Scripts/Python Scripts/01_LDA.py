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

Solution: w = S_W^{-1} (mu_1 - mu_2)   [binary case]

Outputs
-------
lda_results.png      – ROC curve, score distributions, confusion matrix,
                       top-20 gene coefficients
lda_coefficients.csv – full discriminant coefficients for all 81 genes
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
warnings.filterwarnings("ignore")

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.model_selection import StratifiedKFold, cross_val_predict
from sklearn.pipeline import Pipeline
from sklearn.metrics import (
    roc_curve, auc, classification_report,
    confusion_matrix, ConfusionMatrixDisplay,
)

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_PATH = (
    r"C:\Users\Brayan Gutierrez\Desktop\RNAseq-AMD\Dataset\aak100_cpmdat.csv"
)
OUT_DIR = r"C:\Users\Brayan Gutierrez\Desktop\RNAseq-AMD\Scripts\Python Scripts"

# ── 1. Load Data ──────────────────────────────────────────────────────────────
df = pd.read_csv(DATA_PATH, index_col=0)

feature_cols = [c for c in df.columns if c != "mgs_level"]
X = df[feature_cols].values
y = df["mgs_level"].values

le = LabelEncoder()
y_enc = le.fit_transform(y)          # MGS1 → 0,  MGS4 → 1

print("=" * 60)
print("Linear Discriminant Analysis – AMD (MGS1 vs MGS4)")
print("=" * 60)
print(f"Classes  : {le.classes_}")
print(f"Samples  : {len(y_enc)}  (MGS1={( y_enc==0).sum()}, MGS4={(y_enc==1).sum()})")
print(f"Features : {X.shape[1]}")

# ── 2. LDA Pipeline with 10-fold Stratified CV ────────────────────────────────
pipe = Pipeline([
    ("scaler", StandardScaler()),
    ("lda",    LinearDiscriminantAnalysis(solver="svd")),
])

cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

# Cross-validated probability scores (for ROC / metrics)
y_prob = cross_val_predict(pipe, X, y_enc, cv=cv, method="predict_proba")[:, 1]
y_pred = (y_prob >= 0.5).astype(int)

# ── 3. ROC & Metrics ──────────────────────────────────────────────────────────
fpr, tpr, _ = roc_curve(y_enc, y_prob)
roc_auc = auc(fpr, tpr)

print(f"\nCross-validated ROC-AUC : {roc_auc:.4f}")
print("\nClassification Report (10-fold CV):")
print(classification_report(y_enc, y_pred, target_names=le.classes_))

# ── 4. Fit on Full Dataset (for coefficients & projection) ────────────────────
scaler_full = StandardScaler()
X_scaled    = scaler_full.fit_transform(X)

lda_full = LinearDiscriminantAnalysis(solver="svd")
lda_full.fit(X_scaled, y_enc)
X_lda = lda_full.transform(X_scaled)          # (n, 1) projected scores

# Discriminant coefficients (standardised → directly comparable)
coef_df = pd.DataFrame({
    "Gene":        feature_cols,
    "Coefficient": lda_full.coef_[0],
}).sort_values("Coefficient", key=abs, ascending=False).reset_index(drop=True)

print("\nTop 20 discriminant genes (|coefficient|):")
print(coef_df.head(20).to_string(index=False))

# ── 5. Four-panel Figure ─────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(13, 10))
fig.suptitle(
    "Linear Discriminant Analysis – AMD (MGS1 vs MGS4)",
    fontsize=14, fontweight="bold",
)

# Panel A – ROC curve
ax = axes[0, 0]
ax.plot(fpr, tpr, color="steelblue", lw=2.5,
        label=f"LDA  (AUC = {roc_auc:.3f})")
ax.fill_between(fpr, tpr, alpha=0.08, color="steelblue")
ax.plot([0, 1], [0, 1], "k--", lw=1)
ax.set_xlabel("False Positive Rate")
ax.set_ylabel("True Positive Rate")
ax.set_title("(A)  ROC Curve  [10-fold CV]")
ax.legend(loc="lower right", fontsize=10)
ax.grid(alpha=0.3)

# Panel B – Discriminant score distributions
ax = axes[0, 1]
for cls_idx, cls_name in enumerate(le.classes_):
    scores = X_lda[y_enc == cls_idx, 0]
    ax.hist(scores, bins=20, alpha=0.55, density=True,
            label=cls_name, edgecolor="white", linewidth=0.4)
ax.set_xlabel("LD1 Score")
ax.set_ylabel("Density")
ax.set_title("(B)  Discriminant Score Distributions")
ax.legend()
ax.grid(alpha=0.3)

# Panel C – Confusion matrix
ax = axes[1, 0]
cm = confusion_matrix(y_enc, y_pred)
disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=le.classes_)
disp.plot(ax=ax, colorbar=False, cmap="Blues")
ax.set_title("(C)  Confusion Matrix  [10-fold CV]")

# Panel D – Top-20 gene coefficients
ax = axes[1, 1]
top20 = coef_df.head(20).copy()
colors = ["#d62728" if v > 0 else "#1f77b4" for v in top20["Coefficient"]]
ax.barh(
    top20["Gene"][::-1],
    top20["Coefficient"][::-1],
    color=colors[::-1],
    edgecolor="white",
    linewidth=0.4,
)
ax.axvline(0, color="black", linewidth=0.8)
ax.set_xlabel("Discriminant Coefficient (standardised)")
ax.set_title("(D)  Top 20 Genes by |Coefficient|\n"
             "(red = higher in MGS4, blue = higher in MGS1)")
ax.grid(axis="x", alpha=0.3)

plt.tight_layout()
fig.savefig(f"{OUT_DIR}\\lda_results.png", dpi=150, bbox_inches="tight")
plt.show()
print(f"\nFigure saved → lda_results.png")

# ── 6. Save Coefficients ──────────────────────────────────────────────────────
coef_df.to_csv(f"{OUT_DIR}\\lda_coefficients.csv", index=False)
print(f"Coefficients saved → lda_coefficients.csv")
