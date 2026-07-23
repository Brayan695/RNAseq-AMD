# Evaluates the latent embeddings produced by run_cncvae.py: fits simple
# classifiers on each fold's train embedding and scores them on the held-out
# test embedding, so you can see whether the learned representation actually
# separates the MGS classes. Reads the .npz files run_cncvae.py wrote - doesn't
# retrain or touch the VAE itself.
import argparse
import os
import warnings

warnings.filterwarnings('ignore')

import matplotlib

matplotlib.use('Agg')  # headless backend - this script may run without a display
import matplotlib.pyplot as plt
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.manifold import TSNE
from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC

# Classifiers selectable via --classifiers; each is a zero-arg factory so a
# fresh, unfitted instance is created per fold (fitting one classifier on
# multiple folds without resetting it would leak information across folds).
CLASSIFIERS = {
    'logreg': lambda: LogisticRegression(max_iter=1000),
    'nb': lambda: GaussianNB(),
    'svm': lambda: SVC(C=1.5, kernel='rbf', probability=True, random_state=42, gamma='auto'),
    'rf': lambda: RandomForestClassifier(n_estimators=200, random_state=42),
}


def evaluate_fold(emb_train, y_train, emb_test, y_test, classifier_names):
    """Fit each requested classifier on one fold's train embedding, score it on
    that fold's train AND test embedding. Returns {classifier_name: (train_acc,
    test_acc, test_auc)}. AUC assumes a binary label (index [:, 1] = positive
    class probability) - falls back to NaN if that doesn't apply."""
    results = {}
    for name in classifier_names:
        clf = CLASSIFIERS[name]()
        clf.fit(emb_train, y_train)
        train_acc = accuracy_score(y_train, clf.predict(emb_train))
        test_acc = accuracy_score(y_test, clf.predict(emb_test))
        try:
            test_auc = roc_auc_score(y_test, clf.predict_proba(emb_test)[:, 1])
        except Exception:
            test_auc = float('nan')
        results[name] = (train_acc, test_acc, test_auc)
    return results


def main():
    parser = argparse.ArgumentParser(
        description='Evaluate CNC-VAE latent embeddings against the MGS AMD severity label '
                    '(MGS1 vs MGS4), using embeddings produced by run_cncvae.py.'
    )
    parser.add_argument('--resdir', required=True,
                         help='Directory containing the saved embeddings, e.g. '
                              'results/CNCVAE_Clin+mRNA_integration/cncvae_LS_16_DS_128_kl_beta_1.0')
    parser.add_argument('--label_col', default='mgs_level')
    parser.add_argument('--numfolds', type=int, default=5,
                         help='Number of CV folds to aggregate, or 0 to inspect the whole-cohort embedding')
    parser.add_argument('--classifiers', nargs='+', default=['logreg', 'rf'], choices=list(CLASSIFIERS))
    parser.add_argument('--tsne', action='store_true', help='Save a t-SNE plot of the whole-cohort embedding')
    parser.add_argument('--out', default='analyses')
    args = parser.parse_args()

    os.makedirs(args.out, exist_ok=True)

    if args.numfolds == 0:
        # Whole-cohort mode: just load and inspect the single embedding run_cncvae.py
        # --fold 0 produced (no train/test split to score classifiers against).
        embed = np.load(os.path.join(args.resdir, args.label_col + '.npz'))
        labels = np.load(os.path.join(args.resdir, args.label_col + '_labels.npz'), allow_pickle=True)
        emb_train = embed['emb_train']
        y = labels['y']
        classes = labels['classes']

        print('Whole-cohort embedding shape:', emb_train.shape)

        if args.tsne:
            perplexity = max(5, min(30, len(y) - 1))
            proj = TSNE(random_state=42, perplexity=perplexity).fit_transform(emb_train)
            fig, ax = plt.subplots(figsize=(6, 6))
            for cls_idx, cls_name in enumerate(classes):
                mask = y == cls_idx
                ax.scatter(proj[mask, 0], proj[mask, 1], label=cls_name)
            ax.legend()
            ax.set_xlabel('t-SNE 1')
            ax.set_ylabel('t-SNE 2')
            ax.set_title('CNC-VAE latent space')
            fig.savefig(os.path.join(args.out, 'tsne_whole.png'), bbox_inches='tight')
            print('Saved', os.path.join(args.out, 'tsne_whole.png'))
        return

    # Fold mode: load each fold's embedding (saved by run_cncvae.py --fold 1..N)
    # and aggregate classifier performance across all folds into mean +/- std.
    all_results = {name: {'train_acc': [], 'test_acc': [], 'test_auc': []} for name in args.classifiers}
    for fold in range(1, args.numfolds + 1):
        embed = np.load(os.path.join(args.resdir, args.label_col + str(fold) + '.npz'))
        labels = np.load(os.path.join(args.resdir, args.label_col + str(fold) + '_labels.npz'))
        fold_results = evaluate_fold(
            embed['emb_train'], labels['y_train'], embed['emb_test'], labels['y_test'], args.classifiers
        )
        for name, (train_acc, test_acc, test_auc) in fold_results.items():
            all_results[name]['train_acc'].append(train_acc)
            all_results[name]['test_acc'].append(test_acc)
            all_results[name]['test_auc'].append(test_auc)

    out_file = os.path.join(args.out, 'cncvae_{}_results.csv'.format(args.label_col))
    with open(out_file, 'w') as f:
        f.write('classifier,train_acc_mean,train_acc_std,test_acc_mean,test_acc_std,test_auc_mean,test_auc_std\n')
        for name, res in all_results.items():
            f.write('{},{},{},{},{},{},{}\n'.format(
                name,
                np.mean(res['train_acc']), np.std(res['train_acc']),
                np.mean(res['test_acc']), np.std(res['test_acc']),
                np.nanmean(res['test_auc']), np.nanstd(res['test_auc']),
            ))
    print('Wrote', out_file)


if __name__ == '__main__':
    main()
