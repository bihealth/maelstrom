import gc
import multiprocessing
import random
import sqlite3
import sys
import typing

import warnings

from logzero import logger

import attr
from guppy import hpy
import pandas as pd
import numpy as np
import scipy.stats as ss
import tqdm
from sklearn import mixture
from sklearn import ensemble
from sklearn import preprocessing
from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import roc_curve


@attr.s(frozen=True, auto_attribs=True)
class SvFeatureRandomForestConfig:
    """Configuration for ``SvFeatureRandomForest``."""

    max_train_size: int = 100_000
    clean_cutoffs: bool = False

    rf_n_estimators: int = 500
    rf_random_seed: int = 42
    rf_oob_score: bool = True
    rf_max_features: typing.Optional[int] = None


@attr.s(frozen=True, auto_attribs=True)
class _TrainingData:
    """Internal data structure, to hold training data."""

    #: The feature vectors.
    features: pd.DataFrame
    #: Raw labels (NOT encoded yet).
    labels: pd.Series


def _best_cutoff(metric: pd.Series, probs: pd.Series):
    """Compute best cutoff based on a feature and predicted probability vector."""
    preds = metric.values

    # Pass if predicted probability is above 0.5.
    classify = np.vectorize(lambda x: 1 if x >= 0.5 else 0)
    ground_truth = classify(probs)

    # If all variants which passed prior cutoffs also passed random forest,
    # return minimum value instead of trying to compute cutoff
    if 0 not in ground_truth:
        return preds.min()

    fpr, tpr, thresh = roc_curve(ground_truth, preds)
    dist = np.sqrt((fpr - 0) ** 2 + (tpr - 1) ** 2)
    best_idx = np.argmin(dist)

    # If cutoff set at no instances, scikit-learn sets thresh[0] to max(y_score) + 1
    if best_idx == 0:
        return thresh[best_idx] - 1
    else:
        return thresh[best_idx]


class SvFeatureRandomForest:
    """Helper class for running RoC-curve corrected random forest training."""

    def __init__(
        self,
        *,
        df: pd.DataFrame,
        cols_features: typing.Iterable[str],
        col_label: str,
        cols_cutoff_indeps: typing.Iterable[str],
        cols_cutoff_deps: typing.Iterable[str],
        config: typing.Optional[SvFeatureRandomForestConfig] = None,
        label_encoder: typing.Any = None,
    ):
        """
            :df: data frame containing training and prediction data
            :cols_features: names of the feature columns
            :col_label: name of the label column
            :config: configuration of the prediction/scoring
        """
        self.df = df
        self.cols_features = list(cols_features)
        self.col_label = col_label
        self.cols_cutoff_indeps = cols_cutoff_indeps
        self.cols_cutoff_deps = cols_cutoff_deps
        self.config = config or SvFeatureRandomForestConfig()
        self.label_encoder = label_encoder or preprocessing.LabelEncoder().fit([True, False])
        #: Resulting cutoffs, set after calling ``run``.
        self.cutoffs: typing.Optional[typing.Dict[str, float]] = None
        #: Resulting probabilities after applying cutoffs, set after calling ``run``.
        self.probabilites: typing.Optional[pd.Series] = None

    def run(self):
        random_state = np.random.RandomState(self.config.rf_random_seed)
        training_data = self._select_training_data(random_state)
        random_forest = self._fit_random_forest(training_data, random_state)
        test_data = self._select_test_data()
        probs = self._predict_probs(test_data, random_forest)
        if self.config.clean_cutoffs:
            cutoff_data = training_data
        else:
            cutoff_data = test_data
        cutoffs = self._learn_cutoffs(probs, test_data, cutoff_data)
        cutoff_probs = self._cutoff_probs(probs, test_data, cutoffs)
        # Store resulting cutoffs and probabilities.
        self.cutoffs = cutoffs
        self.probabilites = cutoff_probs
        return self

    def _cutoff_probs(self, probs, test_data, cutoffs):
        """Apply cutoffs to probabilities and success/failure enforcement."""
        probs = probs.copy()
        passes_prob_gt_half = pd.Series(probs >= 0.5, name="passes_prob_gt_half")
        passes_all_cutoffs = pd.Series(True, index=probs.index, name="passes_all_cutoffs")
        # Force failure if any metric is below the observed cutoff.
        for metric, value in cutoffs.items():
            passes_cutoff = test_data[metric] >= value
            probs[passes_prob_gt_half & ~passes_cutoff] = 0.499
            passes_all_cutoffs = passes_all_cutoffs & passes_cutoff
        # Forces sucess if all metrics are above the observed cutoff but fails `probs >=0.5`.
        probs[passes_all_cutoffs & ~passes_prob_gt_half] = 0.501
        return probs

    def _learn_cutoffs(
        self, probs: pd.Series, test_data: pd.DataFrame, cutoff_data: pd.DataFrame,
    ) -> typing.Dict[str, float]:
        """Learn the cutoffs based on the data in features."""
        logger.info("Learning cutoffs ...")
        cutoffs = {}
        passing_cutoff = pd.Series(True, index=cutoff_data.index)

        logger.info("  - independent metrics: %s", self.cols_cutoff_indeps)
        for col_name in self.cols_cutoff_indeps:
            feature = cutoff_data[col_name]
            idx = np.searchsorted(test_data.index, cutoff_data.index)
            cutoff = _best_cutoff(feature, probs[idx])
            cutoffs[col_name] = cutoff
            passing_cutoff = passing_cutoff & (cutoff_data[col_name] >= cutoff)

        logger.info("  - dependent metrics: %s", self.cols_cutoff_deps)
        features_passing = cutoff_data.loc[passing_cutoff]
        for col_name in self.cols_cutoff_deps:
            feature = features_passing[col_name]
            idx = np.searchsorted(test_data.index, features_passing.index)
            cutoffs[col_name] = _best_cutoff(feature, probs[idx])

        logger.info("... done learning cutoffs")
        return cutoffs

    def _predict_probs(
        self, test_data: pd.DataFrame, random_forest: ensemble.RandomForestClassifier,
    ) -> pd.Series:
        """Predict probabilities using random forest."""
        logger.info("Predicting probabilities ...")
        probs = random_forest.predict_proba(test_data[self.cols_features].values)
        result = pd.Series(probs[:, 1], index=test_data.index, name="probs")
        logger.info("... done predicting probabilities")
        return result

    def _select_test_data(self) -> pd.DataFrame:
        """Compute tests data from ``self.df``."""
        logger.info("Build test data ...")
        idx_testable = ~self.df[self.cols_features].isnull().any(axis=1)
        result = self.df.loc[idx_testable]
        logger.info("... done building test data")
        return result

    def _fit_random_forest(
        self, training_data: _TrainingData, random_state: np.random.RandomState,
    ) -> pd.DataFrame:
        """Fit the random tree on the training data."""
        logger.info("Fitting random forest ...")
        rf = ensemble.RandomForestClassifier(
            n_estimators=self.config.rf_n_estimators,
            random_state=random_state,
            oob_score=self.config.rf_oob_score,
            max_features=self.config.rf_max_features,
        )
        rf.fit(
            training_data.features.values, self.label_encoder.transform(training_data.labels),
        )
        logger.info("... done fitting random forest")
        return rf

    def _select_training_data(self, random_state: np.random.RandomState,) -> _TrainingData:
        """Build training data with encoded labels."""
        logger.info("Creating training data...")

        # Trainable data must not have any null label.
        idx_labels = ~self.df[self.col_label].isnull()
        idx_features = ~self.df[self.cols_features].isnull().any(axis=1)
        logger.info(
            "  - with cleaned features: %d -> %d", self.df.shape[0], idx_features.shape[0],
        )
        idx_trainable = idx_labels & idx_features
        # Create trainable data with ``bool`` label column.
        df_trainable = self.df.loc[idx_trainable].copy()
        df_trainable[self.col_label] = df_trainable[self.col_label].astype(bool)

        # Downsample the trainable data if necessary.
        if df_trainable.shape[0] >= self.config.max_train_size:
            max_n = self.config.max_train_size // 2
            df_positives = df_trainable[df_trainable[self.col_label]]
            n_pos_orig = df_positives.shape[0]
            if df_positives.shape[0] >= max_n:
                df_positives = df_positives.sample(n=max_n, random_state=random_state)
            df_negatives = df_trainable[~df_trainable[self.col_label]]
            n_neg_orig = df_negatives.shape[0]
            if df_negatives.shape[0] >= max_n:
                df_negatives = df_negatives.sample(n=max_n, random_state=random_state)
            df_trainable = pd.concat(df_positives, df_negatives)
            logger.info(
                "  - downsampled positives %d -> %d, negatives %d -> %d",
                n_pos_orig,
                df_positives.shape[0],
                n_neg_orig,
                df_negatives.shape[0],
            )

        # Perform sanity checks.
        num_pos = df_trainable.loc[df_trainable[self.col_label]].shape[0]
        num_neg = df_trainable.loc[~df_trainable[self.col_label]].shape[0]
        logger.info("  - total positives: %d, negatives: %s", num_pos, num_neg)
        if not num_pos:
            raise Exception("No positive training data in training set")
        if not num_neg:
            raise Exception("No negative training data in training set")

        result = _TrainingData(
            features=df_trainable[self.cols_features], labels=df_trainable[self.col_label],
        )

        logger.info("... done creating training data")
        return result
