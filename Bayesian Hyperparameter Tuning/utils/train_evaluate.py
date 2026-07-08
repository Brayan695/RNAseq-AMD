# This file contains APIs for evaluating models and plotting results

import os
import pickle
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_curve, roc_auc_score, f1_score, precision_recall_curve
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

from xgboost import XGBClassifier
from sklearn.ensemble import RandomForestClassifier

from bayes_opt import BayesianOptimization
from bayes_opt import UtilityFunction

import shap

# Initialize rng
rng = np.random.default_rng(2022)

# apply threshold to positive probabilities to create labels
def to_labels(pos_probs, threshold):
    return (pos_probs >= threshold).astype('int')

def summarize_diagnostics(history):
    '''
    Plots summary diagnostics including loss curves and accuracy vs epochs for any keras model
    
    Inputs:
    history = keras history file generated while training
    
    No outputs
    '''
    fig, ax = plt.subplots(1,2, figsize=(20, 10))
    # plot loss
    ax[0].set_title('Loss Curves', fontsize=20)
    ax[0].plot(history.history['loss'], label='train')
    ax[0].plot(history.history['val_loss'], label='val')
    ax[0].set_xlabel('Epochs', fontsize=15)
    ax[0].set_ylabel('Loss', fontsize=15)
    ax[0].legend(fontsize=15)
    # plot accuracy
    ax[1].set_title('Classification Accuracy', fontsize=20)
    ax[1].plot(history.history['acc'], label='train')
    ax[1].plot(history.history['val_acc'], label='val')
    ax[1].set_xlabel('Epochs', fontsize=15)
    ax[1].set_ylabel('Accuracy', fontsize=15)
    ax[1].legend(fontsize=15)
    
def bayesian_optimization(X,
                          y,
                          n_bootstrap=100,
                          model_type='xgboost',
                          n_iter=100,
                          verbose=1,
                          plot_roc=True,
                          save_roc_fig=True,
                          plot_boxplot=True,
                          save_boxplot=True,
                          save_results=True,
                          dataset_filename='aak81Dataset_cpmdat.csv',
                          path='results'):
    '''
    Runs Bayesian optimization to find the best hyperparameter combination for xgboost and random forest
    
    Inputs:
    X = input features
    y = prediction output
    n_bootstrap = # times to run bootstrap (split data randomly into train, validate, test)
    model_type = 'xgboost' or 'random_forest'
    n_iter = # iterations for Bayesian optimization
    verbose = whether to print iteration # for Bayesian optimization
    plot_roc = whether to plot ROC
    save_roc_fig = whether to save ROC (works only if plot_roc==True)
    plot_boxplot = whether to plot boxplot of sensitivity, specificity, AUROC, precision, recall, and F-1 score
    save_boxplot = wether to save the boxplot
    save_results = whether to save results (saved in a pkl file)
    path = path where figure will be saved
    
    Returns:
    results_df = a pandas data frame with AUROC, Sensitivity, Specificity, F-1 score, Precision, and Recall for n_bootstrap iterations
    best_model_yet = the model with the best test AUROC
    '''
    
    # Build an FPR grid to generate TPR on a grid
    fpr_grid_vec = np.linspace(0, 1, 101, endpoint=True)
    tpr_grid_mat = np.zeros((len(fpr_grid_vec), n_bootstrap))
    
    # Maintain a list of AUROC scores, sensitivity, specificity, precision, recall, and f-1
    auroc_list = []
    sens_list = []
    spec_list = []
    prec_list = []
    rec_list = []
    f1_list = []
    
    # Keep track of best AUROC yet and best model yet
    best_auroc_yet = 0.0
    best_model_yet = None
    
    # Create dataframe of shap values
    shap_df = pd.DataFrame(columns=X.columns)
    
    for i in range(n_bootstrap):
        if verbose:
            print('Bootstrap iteration # {} out of {}'.format(i+1, n_bootstrap))
        # Split into train, validation, and test
        X_train_val, X_test, y_train_val, y_test = train_test_split(X,
                                                                    y,
                                                                    stratify=y,
                                                                    test_size=0.2,
                                                                    random_state=rng.integers(500000))
        X_train, X_val, y_train, y_val = train_test_split(X_train_val,
                                                          y_train_val,
                                                          stratify=y_train_val,
                                                          test_size=0.2,
                                                          random_state=rng.integers(500000))

        if model_type == 'xgboost':
            # Parameter bounds
            param_bounds = {'eta': (0.01, 1.0),
                            'n_estimators': (100, 1000),
                            'gamma':(0.01, 5),
                            'max_depth': (3, 11),
                            'min_child_weight': (1, 11),
                            'subsample': (0.5, 1),
                            'reg_lambda': (0.1, 1000),
                            'reg_alpha': (0.001, 1000),
                            'max_delta_step':(0, 5),
                            'scale_pos_weight':(1.2, 5),
                            'colsample_by': (0.5, 1)}

            # Function for training xgboost given hyperparameters
            def xgboost_hyperparam(eta,
                                   n_estimators,
                                   gamma,
                                   max_depth,
                                   min_child_weight,
                                   subsample,
                                   reg_lambda,
                                   reg_alpha,
                                   max_delta_step,
                                   scale_pos_weight,
                                   colsample_by):

                max_depth = int(max_depth)
                n_estimators = int(n_estimators)
                min_child_weight = int(min_child_weight)

                model = XGBClassifier(eta=eta,
                                      n_estimators=n_estimators,
                                      gamma=gamma,
                                      max_depth=max_depth,
                                      min_child_weight=min_child_weight,
                                      subsample=subsample,
                                      reg_lambda=reg_lambda,
                                      reg_alpha=reg_alpha,
                                      max_delta_step=max_delta_step,
                                      scale_pos_weight=scale_pos_weight,
                                      colsample_by=colsample_by,
                                      verbosity=0)

                # Fit on training data
                model.fit(X_train, y_train)

                # Get predictions on validation set
                y_pred = model.predict_proba(X_val)

                # Calculate AUROC
                auroc = roc_auc_score(y_val, y_pred[:,1])
                return auroc

            # Initialize optimizer
            optimizer = BayesianOptimization(f=xgboost_hyperparam,
                                             pbounds=param_bounds,
                                             random_state=1,
                                             verbose=0)

        elif model_type == 'random_forest':
            # Parameter bounds
            param_bounds = {'n_estimators': (100, 1000),
                            'max_depth': (5, 101),
                            'min_samples_split': (2, 11),
                            'min_samples_leaf': (2, 11),
                            'max_features_choose': (0, 2)}
            max_features_list = ['sqrt', 'log2']
            # Function for training random forest given hyperparameters
            def randomForest_hyperparam(n_estimators,
                                        max_depth,
                                        min_samples_split,
                                        min_samples_leaf,
                                        max_features_choose):

                n_estimators = int(n_estimators)
                max_depth = int(max_depth)
                min_samples_split = int(min_samples_split)
                min_samples_leaf = int(min_samples_leaf)
                max_features = max_features_list[int(max_features_choose-0.0000001)]


                model = RandomForestClassifier(n_estimators=n_estimators,
                                               max_depth=max_depth,
                                               min_samples_split=min_samples_split,
                                               min_samples_leaf=min_samples_leaf,
                                               max_features=max_features,
                                               n_jobs=1000,
                                               random_state=rng.integers(500000))

                # Fit on training data
                model.fit(X_train, y_train)

                # Get predictions on validation set
                y_pred = model.predict_proba(X_val)

                # Calculate AUROC
                auroc = roc_auc_score(y_val, y_pred[:,1])
                return auroc

            # Initialize optimizer
            optimizer = BayesianOptimization(f=randomForest_hyperparam,
                                             pbounds=param_bounds,
                                             random_state=1,
                                             verbose=0)

        # Set Gaussian Process parameters here
        optimizer.set_gp_params(alpha=1e-5)

        # Define acquisition strategy here
        utility = UtilityFunction(kind="poi", xi=0.0)
        
        # Run optimizer
        optimizer.maximize(init_points=20,
                           n_iter=n_iter,
                           acquisition_function=utility)
        
        # Get best hyperparameter combination from optimizer results
        best_auroc = 0.
        best_params = None
        for details in optimizer.res:
            auroc = details['target']
            if auroc > best_auroc:
                best_auroc = auroc
                best_params = details
                
        # Create model with best hyperparams
        if model_type == 'xgboost':
            best_xgb_params = dict(best_params['params'])
            best_xgb_params['max_depth'] = int(best_xgb_params['max_depth'])
            best_xgb_params['n_estimators'] = int(best_xgb_params['n_estimators'])
            best_xgb_params['min_child_weight'] = int(best_xgb_params['min_child_weight'])
            model = XGBClassifier(**best_xgb_params, verbosity=0)
        elif model_type == 'random_forest':
            model = RandomForestClassifier(n_estimators=int(best_params['params']['n_estimators']),
                                           max_depth=int(best_params['params']['max_depth']),
                                           min_samples_split=int(best_params['params']['min_samples_split']),
                                           min_samples_leaf=int(best_params['params']['min_samples_leaf']),
                                           max_features=max_features_list[int(best_params['params']['max_features_choose']-0.0000001)])
            
        # Fit on training data
        model.fit(X_train, y_train)
        
        # Get SHAP values
        explainer = shap.TreeExplainer(model)
        shap_values = explainer.shap_values(X)
        
        # Append SHAP values to dataframe
        shap_df = pd.concat([shap_df, pd.DataFrame(shap_values, columns=shap_df.columns)], ignore_index=True)

        # Predict
        y_pred = model.predict_proba(X_test)

        fpr, tpr, thresholds = roc_curve(y_true=y_test, y_score=y_pred[:,1])

        # get the best threshold in ROC using Youden's J statistic
        J = tpr - fpr
        ix = np.argmax(J)
        best_thresh = thresholds[ix]

        # Get balanced values of fpr and tpr
        best_fpr = fpr[ix]
        best_tpr = tpr[ix]

        # Get best balanced values of sensitivity and specificity
        best_sens = best_tpr
        best_spec = 1 - best_fpr

        # Append to list
        sens_list.append(best_sens)
        spec_list.append(best_spec)

        # Get best threshold for balance of precision and recall using F-1 score
        precision, recall, thresholds = precision_recall_curve(y_test, y_pred[:,1])
        # convert to f score
        f1score = (2 * precision * recall) / (precision + recall)
        # locate the index of the largest f score
        ix = np.argmax(f1score)
        # Get best balanced precision and recall using F-1 score
        best_precision = precision[ix]
        best_recall = recall[ix]
        # Append to list
        prec_list.append(best_precision)
        rec_list.append(best_recall)
        f1_list.append(f1score[ix])

        # Interpolation for plotting
        fpr = np.concatenate(([0], fpr, [1]))
        tpr = np.concatenate(([0], tpr, [1]))

        interpolator = interp1d(fpr, tpr, kind='nearest')

        tpr_grid_mat[:, i] = interpolator(fpr_grid_vec)

        auroc = roc_auc_score(y_test, y_pred[:,1])
        auroc_list.append(auroc)

        if auroc > best_auroc_yet:
            best_auroc_yet = auroc
            best_model_yet = model
            
    # confidence interval for AUCs
    sorted_scores = np.array(auroc_list)
    sorted_scores.sort()
    auroc_mean = np.mean(sorted_scores)
    confidence_lower = sorted_scores[int(0.05 * len(sorted_scores))]
    confidence_upper = sorted_scores[int(0.95 * len(sorted_scores))]
    print('AUROC: {} [CI {} - {}]'.format(auroc_mean, confidence_lower, confidence_upper))
    
    # Prefix for saving if needed
    if model_type == 'xgboost':
        prefix = 'XGBoost'
    elif model_type == 'random_forest':
        prefix = 'RandomForest'

    if plot_roc:
        # Plot
        # confidence interval for ROC
        tpr_grid_mats = np.sort(tpr_grid_mat, axis=1)
        tpr_low_05 = tpr_grid_mats[:, int(0.05 * n_bootstrap)]
        tpr_top_95 = tpr_grid_mats[:, int(0.95 * n_bootstrap)]
        tpr_mean = np.mean(tpr_grid_mat, axis=1)

        plt.figure(figsize=(10,10))
        plt.xlabel('True Negative Rate', fontsize=24)
        plt.ylabel('True Positive Rate', fontsize=24)
        plt.tick_params(axis='both', which='major', labelsize=18)
        ax = plt.gca()  # kwargs.pop('ax', plt.gca())
        
        base_line, = ax.plot(fpr_grid_vec,
                             tpr_mean,
                             '-',
                             linewidth=4,
                             label=prefix+'\nAUC: {:.2f} [{:.2f} - {:.2f}]'.format(auroc_mean,
                                                                                   confidence_lower,
                                                                                   confidence_upper))
        ax.fill_between(fpr_grid_vec, tpr_low_05, tpr_top_95, facecolor=base_line.get_color(), alpha=0.2);
        plt.legend(loc='lower right', fontsize=18);
        
        if save_roc_fig:
            plt.savefig(os.path.join(path, '{}_Bayesian_ROC_{}.pdf'.format(prefix, dataset_filename)))
            
    if save_results:
        # Save results in pkl file
        results_dict = {'auroc_mean': auroc_mean,
                        'confidence_lower': confidence_lower,
                        'confidence_upper': confidence_upper,
                        'auroc_list': auroc_list,
                        'fpr_grid_vec': fpr_grid_vec,
                        'tpr_grid_mats': tpr_grid_mats,
                        'tpr_low_05': tpr_low_05,
                        'tpr_top_95': tpr_top_95,
                        'tpr_mean': tpr_mean,
                        'sens_list': sens_list,
                        'spec_list': spec_list,
                        'prec_list': prec_list,
                        'rec_list': rec_list,
                        'f1_list': f1_list,
                        }
        with open(os.path.join(path, '{}_Bayesian_results_{}.pkl'.format(prefix,
                                                                         dataset_filename)), 'wb') as f:
            pickle.dump(results_dict, f)
            
        shap_df.to_csv(os.path.join(path, '{}_Bayesian_results_shap_{}.csv'.format(prefix,
                                                                                   dataset_filename)))
            
            
    # Prepare results to return
    results_df = pd.DataFrame(list(zip(auroc_list,
                                   sens_list,
                                   spec_list,
                                   f1_list,
                                   prec_list,
                                   rec_list)),
                         columns=['AUROC',
                                  'Sensitivity',
                                  'Specificity',
                                  'F-1 score',
                                  'Precision',
                                  'Recall'])
    
    if plot_boxplot:
        fig = plt.figure(figsize=(12,12))
        ax = fig.add_subplot(1,1,1)
        sns.boxplot(data=results_df, ax=ax)
        plt.tick_params(axis='both', which='major', labelsize=18)
        
        # Save if needed
        if save_boxplot:
            plt.savefig(os.path.join(path, '{}_Bayesian_boxplot_{}.pdf'.format(prefix,
                                                                               dataset_filename)))
    
    return results_df, best_model_yet, shap_df

    
def mlp_bootstrap(n_bootstrap = 10,
                  n_iter = 10):
    '''
    Trains and evaluates MLP model with bootstrap, Plots ROC with CIs, Prints AUROC with CI
    
    Inputs:
    n_bootstraps = how many times to run bootstrap (split data randomly into train, validate, test)
    n_iter = how many times to instantiate and train a model for a given split
    
    
    
    '''
    
    # Build an FPR grid to generate TPR on a grid
    fpr_grid_vec = np.linspace(0, 1, 100, endpoint=True)
    tpr_grid_mat = np.zeros((len(fpr_grid_vec), n_bootstrap*n_iter))

    # Maintain a list of AUROC scores
    auroc_list = []

    # index variable
    count = 0

    # Run experiments
    for i in range(n_bootstrap):
        print(i)
        # Split data into train, validation, and test
        # train = 64%, val = 16%, test = 20%
        X_train_val, X_test, y_train_val, y_test = train_test_split(X,
                                                                    y,
                                                                    stratify=y,
                                                                    test_size=0.2,
                                                                    random_state=rng.integers(500000))
        X_train, X_val, y_train, y_val = train_test_split(X_train_val,
                                                        y_train_val,
                                                        stratify=y_train_val,
                                                        test_size=0.2,
                                                        random_state=rng.integers(500000))

        for j in range(n_iter):
            # define model
            tf.random.set_seed(rng.integers(500000))
            opt = Adam(learning_rate=0.0001)
            input_shape = X_train.shape[1]
            model = mlp_model(input_shape,
                              n_layers=2,
                              n_neurons=[64, 64],
                              opt=opt,
                              dropout=False,
                              dropout_rate=0.2,
                              batchnorm=False)
        #     model.summary()

            # Early stopping
            es = EarlyStopping(monitor='val_loss',
                               mode='min',
                               verbose=0,
                               patience=50,
                               restore_best_weights=True)

            # Train a model
            history = model.fit(X_train,
                                y_train,
                                epochs=200,
                                batch_size=16,
                                validation_data=(X_val, y_val),
                                verbose=0,
                                callbacks=[es])

            # Plot loss curves
        #     summarize_diagnostics(history)

            # Predict on test set
            y_pred = model.predict(X_test)

            # Measure ROC
        #     roc_auc = roc_auc_score(y_test, y_pred)
        #     print('AUROC:',roc_auc)
        #     test_acc = accuracy_score(y_test, np.around(y_pred))
        #     print('accuracy:', test_acc)
        #     p = precision_score(y_test, np.around(y_pred))
        #     print('Precision:', p)
        #     r = recall_score(y_test, np.around(y_pred))
        #     print('Recall:', r)
        #     confusionMatrix = confusion_matrix(y_test, np.around(y_pred))
        #     print('Confusion Matrix:\n', confusionMatrix)
            fpr, tpr, _ = roc_curve(y_true=y_test, y_score=y_pred)
            fpr = np.concatenate(([0], fpr, [1]))
            tpr = np.concatenate(([0], tpr, [1]))

            interpolator = interp1d(fpr, tpr, kind='nearest')

            tpr_grid_mat[:, count] = interpolator(fpr_grid_vec)

            auroc_list.append(roc_auc_score(y_test, y_pred))

            count += 1
    
    # confidence interval for AUCs
    sorted_scores = np.array(auroc_list)
    sorted_scores.sort()
    auroc_mean = np.mean(sorted_scores)
    confidence_lower = sorted_scores[int(0.025 * len(sorted_scores))]
    confidence_upper = sorted_scores[int(0.975 * len(sorted_scores))]
    print('AUROC: {} [CI {} - {}]'.format(auroc_mean, confidence_lower, confidence_upper))

    # Plot
    # confidence interval for ROC
    tpr_grid_mats = np.sort(tpr_grid_mat, axis=1)
    tpr_low_025 = tpr_grid_mats[:, int(0.025 * n_bootstrap * n_iter)]
    tpr_top_975 = tpr_grid_mats[:, int(0.975 * n_bootstrap * n_iter)]
    tpr_mean = np.mean(tpr_grid_mat, axis=1)

    # plt.hold(True)
    ax = plt.gca()  # kwargs.pop('ax', plt.gca())
    base_line, = ax.plot(fpr_grid_vec, tpr_mean, '-', linewidth=4)
    ax.fill_between(fpr_grid_vec, tpr_low_025, tpr_top_975, facecolor=base_line.get_color(), alpha=0.2)