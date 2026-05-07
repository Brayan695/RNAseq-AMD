"""
Kernel Fisher Discriminant Analysis (KFDA) for AMD Classification
==================================================================
Dataset : aak100_cpmdat.csv  (CPM-normalised RNA-seq, 81 genes, 166 samples)
Classes : MGS1 (early/normal, n=105) vs MGS4 (advanced AMD, n=61)

Method reference
----------------
Mika et al. (1999). "Fisher Discriminant Analysis with Kernels."
See also: https://en.wikipedia.org/wiki/Kernel_Fisher_discriminant_analysis

KFDA implicitly maps data into a (possibly infinite-dimensional) feature
space via a kernel function k(x,y) = <phi(x), phi(y)>, then applies
Fisher's criterion in that space.

Between-class scatter (kernel space):
    M = sum_c  n_c (m_c - m)(m_c - m)^T

Within-class scatter (kernel space, with Tikhonov regularisation):
    N = sum_c  K_c^T K_c  +  epsilon * I

Projection:  alpha = N^{-1} (m_1 - m_2)    [binary case]
Score for a test point x:  f(x) = sum_i alpha_i k(x_i, x)

Kernels evaluated
-----------------
  RBF  (gamma = 0.01, 0.001)
  Polynomial (degree = 2)
  Linear (recovers standard Fisher LDA in primal)

Outputs
-------
kfda_results.png – ROC curves for all kernels, score distribution (best
                   kernel), confusion matrix, AUC bar chart
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings("ignore")

from sklearn.base import BaseEstimator, ClassifierMixin
from sklearn.preprocessing import LabelEncoder, StandardScaler
from sklearn.model_selection import StratifiedKFold
from sklearn.metrics import (
    roc_curve, auc, classification_report,
    confusion_matrix, ConfusionMatrixDisplay,
)
from scipy.linalg import eigh

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_PATH = (
    r"C:\Users\Brayan Gutierrez\Desktop\RNAseq-AMD\Dataset\aak100_cpmdat.csv"
)
OUT_DIR = r"C:\Users\Brayan Gutierrez\Desktop\RNAseq-AMD\Scripts\Python Scripts"


# ══════════════════════════════════════════════════════════════════════════════
# KFDA Classifier
# ══════════════════════════════════════════════════════════════════════════════
class KernelFisherDA(BaseEstimator, ClassifierMixin):
    """
    Binary Kernel Fisher Discriminant Analysis.

    Parameters
    ----------
    kernel : {'rbf', 'poly', 'linear'}
    gamma  : float  – RBF bandwidth parameter
    degree : int    – polynomial degree
    coef0  : float  – polynomial bias term
    reg    : float  – Tikhonov regularisation added to within-class scatter N
    """

    def __init__(self, kernel="rbf", gamma=0.01, degree=2,
                 coef0=1.0, reg=1e-2):
        self.kernel = kernel
        self.gamma  = gamma
        self.degree = degree
        self.coef0  = coef0
        self.reg    = reg

    # ── Kernel computation ────────────────────────────────────────────────────
    def _K(self, A, B):
        if self.kernel == "rbf":
            # ||a - b||^2 via outer products (numerically stable)
            sq = (
                np.sum(A ** 2, axis=1, keepdims=True)
                + np.sum(B ** 2, axis=1, keepdims=True).T
                - 2.0 * A @ B.T
            )
            sq = np.maximum(sq, 0.0)          # numerical guard
            return np.exp(-self.gamma * sq)
        elif self.kernel == "poly":
            return (self.coef0 + A @ B.T) ** self.degree
        elif self.kernel == "linear":
            return A @ B.T
        else:
            raise ValueError(f"Unknown kernel: {self.kernel!r}")

    # ── Fit ───────────────────────────────────────────────────────────────────
    def fit(self, X, y):
        self.X_tr_     = X.copy()
        self.classes_  = np.unique(y)
        assert len(self.classes_) == 2, "KFDA supports binary classification."

        n        = X.shape[0]
        K        = self._K(X, X)              # (n, n) raw kernel matrix

        # ── Centre kernel matrix in feature space ─────────────────────────────
        # K_c = K - 1_n K - K 1_n + 1_n K 1_n
        row_mean  = K.mean(axis=1, keepdims=True)   # (n, 1)
        col_mean  = K.mean(axis=0, keepdims=True)   # (1, n)
        grand_mean = K.mean()

        K_c = K - row_mean - col_mean + grand_mean

        # Cache for later centring of test points
        self._col_mean_   = col_mean
        self._grand_mean_ = grand_mean

        # ── Class indices & centred class means ───────────────────────────────
        i0 = np.where(y == self.classes_[0])[0]
        i1 = np.where(y == self.classes_[1])[0]
        n0, n1 = len(i0), len(i1)

        m0  = K_c[i0].mean(axis=0)    # mean kernel-row for class 0  (n,)
        m1  = K_c[i1].mean(axis=0)    # mean kernel-row for class 1  (n,)
        m   = K_c.mean(axis=0)        # overall mean                  (n,)

        # ── Between-class scatter M ───────────────────────────────────────────
        d0 = (m0 - m)[:, None]        # (n, 1)
        d1 = (m1 - m)[:, None]
        M  = n0 * d0 @ d0.T + n1 * d1 @ d1.T   # (n, n)

        # ── Within-class scatter N  (+ regularisation) ───────────────────────
        K0 = K_c[i0] - m0             # (n0, n)  centred class rows
        K1 = K_c[i1] - m1             # (n1, n)
        N  = K0.T @ K0 + K1.T @ K1 + self.reg * np.eye(n)

        # ── Generalised eigenvalue problem:  M alpha = lambda N alpha ─────────
        vals, vecs = eigh(M, N)
        order = np.argsort(vals)[::-1]
        self.alphas_      = vecs[:, order]       # columns = projection directions
        self.eigenvalues_ = vals[order]

        # ── Projection of training data (for decision boundary) ───────────────
        self._K_c_train_ = K_c
        self._i0_, self._i1_ = i0, i1
        return self

    # ── Centre a test kernel matrix ───────────────────────────────────────────
    def _centre_test(self, K_test):
        """K_test : (n_test, n_train) raw kernel matrix."""
        return (
            K_test
            - K_test.mean(axis=1, keepdims=True)
            - self._col_mean_
            + self._grand_mean_
        )

    # ── Project to discriminant subspace ─────────────────────────────────────
    def _project(self, X_new, n_components=1):
        K_test   = self._K(X_new, self.X_tr_)
        K_test_c = self._centre_test(K_test)
        return K_test_c @ self.alphas_[:, :n_components]

    # ── Decision function: signed distance from midpoint of class centroids ───
    def decision_function(self, X):
        Z      = self._project(X, n_components=1).ravel()
        Z_tr   = (self._K_c_train_ @ self.alphas_[:, :1]).ravel()
        c0_val = Z_tr[self._i0_].mean()
        c1_val = Z_tr[self._i1_].mean()
        sign   = np.sign(c1_val - c0_val)   # ensure MGS4 scores higher than MGS1
        return sign * (Z - 0.5 * (c0_val + c1_val))

    # ── Soft probabilities via logistic calibration ───────────────────────────
    def predict_proba(self, X):
        d     = self.decision_function(X)
        prob1 = 1.0 / (1.0 + np.exp(-d))
        return np.column_stack([1.0 - prob1, prob1])

    def predict(self, X):
        return self.classes_[(self.decision_function(X) >= 0).astype(int)]


# ══════════════════════════════════════════════════════════════════════════════
# Main Analysis
# ══════════════════════════════════════════════════════════════════════════════

# ── 1. Load Data ──────────────────────────────────────────────────────────────
df = pd.read_csv(DATA_PATH, index_col=0)

feature_cols = [c for c in df.columns if c != "mgs_level"]
X = df[feature_cols].values
y = df["mgs_level"].values

le    = LabelEncoder()
y_enc = le.fit_transform(y)

print("=" * 60)
print("Kernel Fisher Discriminant Analysis – AMD (MGS1 vs MGS4)")
print("=" * 60)
print(f"Classes  : {le.classes_}")
print(f"Samples  : {len(y_enc)}  (MGS1={(y_enc==0).sum()}, MGS4={(y_enc==1).sum()})")
print(f"Features : {X.shape[1]}")

scaler   = StandardScaler()
X_scaled = scaler.fit_transform(X)

# ── 2. Cross-Validation for Each Kernel ──────────────────────────────────────
kernels = {
    "RBF (γ=0.01)":     KernelFisherDA(kernel="rbf",    gamma=0.01,  reg=1e-2),
    "RBF (γ=0.001)":    KernelFisherDA(kernel="rbf",    gamma=0.001, reg=1e-2),
    "Polynomial (d=2)": KernelFisherDA(kernel="poly",   degree=2,    reg=1e-2),
    "Linear":           KernelFisherDA(kernel="linear",              reg=1e-2),
}

cv      = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)
results = {}

print("\n10-fold CV results:")
for name, model in kernels.items():
    y_probs = np.zeros(len(y_enc))
    y_preds = np.zeros(len(y_enc), dtype=int)

    for tr, te in cv.split(X_scaled, y_enc):
        m = model.__class__(**model.get_params())
        m.fit(X_scaled[tr], y_enc[tr])
        p_hat          = m.predict_proba(X_scaled[te])[:, 1]
        y_probs[te]    = p_hat
        y_preds[te]    = (p_hat >= 0.5).astype(int)

    fpr, tpr, _ = roc_curve(y_enc, y_probs)
    roc_auc     = auc(fpr, tpr)

    results[name] = {
        "fpr": fpr, "tpr": tpr, "auc": roc_auc,
        "y_prob": y_probs, "y_pred": y_preds,
    }
    print(f"  {name:<25s}  AUC = {roc_auc:.4f}")

# ── 3. Best Kernel ────────────────────────────────────────────────────────────
best_name = max(results, key=lambda k: results[k]["auc"])
print(f"\nBest kernel : {best_name}  (AUC = {results[best_name]['auc']:.4f})")
print("\nClassification Report (10-fold CV, best kernel):")
print(classification_report(
    y_enc, results[best_name]["y_pred"], target_names=le.classes_
))

# ── 4. Fit Best Model on Full Data (for score distributions) ─────────────────
best_params = kernels[best_name].get_params()
best_model  = KernelFisherDA(**best_params)
best_model.fit(X_scaled, y_enc)
Z_full = best_model._project(X_scaled, n_components=1).ravel()

# ── 5. Four-panel Figure ──────────────────────────────────────────────────────
PALETTE = ["#1f77b4", "#ff7f0e", "#2ca02c", "#d62728"]

fig, axes = plt.subplots(2, 2, figsize=(13, 10))
fig.suptitle(
    "Kernel Fisher Discriminant Analysis – AMD (MGS1 vs MGS4)",
    fontsize=14, fontweight="bold",
)

# Panel A – ROC curves (all kernels)
ax = axes[0, 0]
for (name, res), col in zip(results.items(), PALETTE):
    ax.plot(res["fpr"], res["tpr"], lw=2, color=col,
            label=f"{name}  (AUC={res['auc']:.3f})")
ax.plot([0, 1], [0, 1], "k--", lw=1)
ax.set_xlabel("False Positive Rate")
ax.set_ylabel("True Positive Rate")
ax.set_title("(A)  ROC Curves  [10-fold CV]")
ax.legend(loc="lower right", fontsize=8)
ax.grid(alpha=0.3)

# Panel B – Score distributions (best kernel)
ax = axes[0, 1]
for cls_idx, cls_name in enumerate(le.classes_):
    ax.hist(Z_full[y_enc == cls_idx], bins=20, alpha=0.55, density=True,
            label=cls_name, edgecolor="white", linewidth=0.4)
ax.set_xlabel("KFD Score")
ax.set_ylabel("Density")
ax.set_title(f"(B)  Score Distributions\n[Best kernel: {best_name}]")
ax.legend()
ax.grid(alpha=0.3)

# Panel C – Confusion matrix (best kernel)
ax = axes[1, 0]
cm   = confusion_matrix(y_enc, results[best_name]["y_pred"])
disp = ConfusionMatrixDisplay(confusion_matrix=cm, display_labels=le.classes_)
disp.plot(ax=ax, colorbar=False, cmap="Blues")
ax.set_title(f"(C)  Confusion Matrix\n[{best_name}, 10-fold CV]")

# Panel D – AUC bar chart
ax = axes[1, 1]
names_list = list(results.keys())
aucs_list  = [results[n]["auc"] for n in names_list]
bars = ax.bar(range(len(names_list)), aucs_list,
              color=PALETTE, alpha=0.82, edgecolor="white")
ax.set_xticks(range(len(names_list)))
ax.set_xticklabels(names_list, rotation=22, ha="right", fontsize=9)
ax.set_ylabel("AUC")
ax.set_ylim(0.45, 1.05)
ax.set_title("(D)  AUC Comparison Across Kernels")
ax.axhline(0.5, color="gray", linestyle="--", linewidth=1)
for bar, val in zip(bars, aucs_list):
    ax.text(bar.get_x() + bar.get_width() / 2,
            bar.get_height() + 0.01,
            f"{val:.3f}", ha="center", va="bottom", fontsize=9)
ax.grid(axis="y", alpha=0.3)

plt.tight_layout()
fig.savefig(f"{OUT_DIR}\\kfda_results.png", dpi=150, bbox_inches="tight")
plt.show()
print(f"\nFigure saved → kfda_results.png")
