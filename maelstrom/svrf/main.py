"""SV classification using random forest.

This is adapted from the ideas by Werling et al. (2018) and Collins et a. (2019).

- Werling, Donna M., et al. "An analytical framework for whole-genome sequence
  association studies and its implications for autism spectrum disorder."
  Nature genetics 50.5 (2018): 727-736.
- Collins, Ryan L., et al. "An open resource of structural variation for
  medical and population genetics." BioRxiv (2019): 578674.
"""

import pathlib
import gc
import multiprocessing as mp
from functools import partial, reduce
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
import statsmodels.stats.power as sp
from sklearn import mixture
from sklearn import ensemble
from sklearn import preprocessing
from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import roc_curve

from . import random_forest

#: Types of SVS
SV_TYPES = ("DEL", "DUP", "INV", "BND")


AUTOSOMES = set(map(str, range(1, 23)))
ALLOSOMES = set(["X", "Y"])


def _fix_df(df):
    """Helper to get rid of NaN/+inf/-inf as 0."""
    df = df.replace([np.inf, -np.inf], np.nan).copy()
    df.fillna(0, inplace=True)
    return df


def run_or_load(name, token, func, args=None, kwargs=None):
    args = args or []
    kwargs = kwargs or {}
    p = pathlib.Path("%s.%s.pkl.gz" % (name, token))
    if p.exists():
        logger.info("Loading %s from %s ...", name, p)
        return pd.read_pickle(str(p))
        logger.info("... done loading %s from %s", name, p)
    else:
        res = func(*args, **kwargs)
        logger.info("Writing %s to %s ...", name, p)
        res.to_pickle(str(p))
        logger.info("... done writing %s to %s", name, p)
        return res


@attr.s(frozen=True, auto_attribs=True)
class Config:
    """Program configuration."""

    #: Number of processors to use.
    num_procs: int = 16
    #: Random seed to use.
    random_seed: int = 42
    #: Number of rows to process from the input.
    # nrows: typing.Optional[int] = None
    nrows: typing.Optional[int] = 10_000

    min_controls: int = 20
    expected_snv_ratio: float = 0.0005

    #: Input pattern.
    input_pattern = (
        "/fast/projects/medgen_genomes/work/2019-06-05_genomes_reboot/"
        "GRCh37/wgs_sv_calling/maelstrom/merged.%s.tsv"
    )
    #: PESR callers
    pesr_callers = ("delly", "manta")
    #: DoC callers
    doc_callers = ("cnmops",)
    #: FORMAT fields of interest
    formats = ("cnmops", "delly", "manta", "BF", "GT", "PR", "RD", "ROH", "SR", "VL", "VM", "VR")
    #: integer typed FORMAT keys
    formats_int = ("cnmops", "delly", "manta", "PR", "ROH", "SR", "VL", "VM", "VR")
    #: float typed FORMAT keys
    formats_float = ("BF", "RD")
    #: FORMATs to normalize by autosome DoC
    doc_normalize_formats = ()  # ("PR", "SR")

    # Path to SVNs BAF sqlite database.
    path_baf_snvs: str = (
        "/fast/work/projects/medgen_genomes/2019-06-05_genomes_reboot/GRCh37/"
        "wgs_sv_calling/maelstrom/work/maelstrom_baf_gather/delly/all.db"
    )

    cnv_min_reliable_size = 5000
    # Humongous CNVs to ignore CNVs for.
    cnv_max_pesr_cnv_len = 10_000_000

    gmm_n_components: int = 3
    gmm_covariance_type: str = "spherical"

    #: Max. number of samples in pairwise Kolmogorov-Smirnov test.
    ks_max_values = 1000

    run_gc: bool = False

    rds_fail_min: float = 0.0
    rds_fail_max: float = 0.15
    rds_pass_min: float = 0.4
    rds_pass_max: float = 1.0

    path_doc_summary: str = (
        "/fast/work/projects/medgen_genomes/2019-06-05_genomes_reboot/GRCh37/"
        "wgs_sv_calling/maelstrom/work/maelstrom_doc_summary/delly/all.tsv"
    )


def load_autosome_doc(config: Config) -> typing.Dict[str, float]:
    logger.info("Loading autosome DoC ...")
    df = pd.read_csv(
        config.path_doc_summary, sep="\t", na_values=("", ".", "nan"), index_col="CHROM",
    )
    columns = [c.replace("-N1-DNA1-WGS1", "") for c in df.columns]
    result = dict(zip(columns, df.loc[:"autosomes"].values[0]))
    logger.info("... done loading autosome DoC")
    return result


@attr.s(frozen=True, auto_attribs=True)
class RawData:
    """Raw data as read from the disk."""

    #: The INFO field values.
    infos: pd.DataFrame
    #: The FORMAT key values by field name.
    formats: typing.Dict[str, pd.DataFrame]
    #: Autosome DoC by sample.
    autosome_doc: typing.Dict[str, float]


@attr.s(frozen=True, auto_attribs=True)
class BafRecord:
    sv_id: str
    sample: str
    bafs: typing.List[float]

    @classmethod
    def construct(cls, **kwargs):
        if isinstance(kwargs["bafs"], list):
            bafs = kwargs["bafs"]
        else:
            bafs = kwargs["bafs"].split(";")
        return BafRecord(
            sample=kwargs["sample"].replace("-N1-DNA1-WGS1", ""), sv_id=kwargs["svid"], bafs=bafs,
        )


def load_sqlite_records(db_path, sv_id):
    def dict_factory(cursor, row):
        d = {}
        for idx, col in enumerate(cursor.description):
            d[col[0]] = row[idx]
        return d

    with sqlite3.connect(db_path) as conn:
        conn.row_factory = dict_factory
        c = conn.cursor()
        c.execute("SELECT * from baf_snvs WHERE svid = ?", (sv_id,))
        for v in c.fetchall():
            if v["bafs"]:
                bafs = list(map(float, v["bafs"].split(";")))
            else:
                bafs = []
            yield BafRecord.construct(
                svid=v["svid"], sample=v["sample"], bafs=bafs,
            )


def load_bafs(
    sv_id: str, sample_pos: typing.Iterable[str], sample_neg: typing.Iterable[str], config: Config
) -> typing.Tuple[typing.List[float], typing.List[float]]:
    """Load BAFs and split by positive and negative samples."""
    baf_records = list(load_sqlite_records(config.path_baf_snvs, sv_id))
    if not baf_records:
        return None, None
    # logger.debug("loaded %d records for %s", len(baf_records), sv_id)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        baf_median = np.median([baf for r in baf_records for baf in r.bafs])
    scaled_records = [
        attr.evolve(r, bafs=[baf / baf_median * 0.5 for baf in r.bafs]) for r in baf_records
    ]

    pos = []
    neg = []
    for record in scaled_records:
        if record.sample in sample_pos:
            pos += record.bafs
        elif record.sample in sample_neg:
            neg += record.bafs

    return pos, neg


def process_duplications(
    raw_data: RawData, config: Config, random_state: np.random.RandomState, caller: str = "delly"
) -> pd.DataFrame:
    logger.info("Selecting data for duplications...")
    # DUPlications.
    idx_dup = raw_data.infos.SVTYPE == "DUP"
    # CNVs above reliable threshold (e.g., 5k) where read-depth signal is reliable.
    idx_cnv_reliable = raw_data.infos.SVTYPE.isin(("DEL", "DUP")) & (
        raw_data.infos.SVLEN >= config.cnv_min_reliable_size
    )
    # CNVs above a certain size are ignored (unreliable by PESR signal).
    idx_cnv_large = raw_data.infos.SVTYPE.isin(("DEL", "DUP")) & (
        raw_data.infos.SVLEN >= config.cnv_max_pesr_cnv_len
    )

    # Select the variants that we consider reliable.
    idx = idx_dup & idx_cnv_reliable & ~idx_cnv_large

    # Variants called by ``caller`` and not called/background.
    mask_calls = (raw_data.formats[caller] != 0)[idx]
    mask_controls = (raw_data.formats[caller] == 0)[idx]
    logger.info("... done selecting data for duplicatios")

    logger.info("Performing pairwise KS tests ...")

    def sample_vals(vals):
        return random_state.choice(vals, min(len(vals), config.ks_max_values), replace=False)

    ks_ids = []
    ks_stats = []
    ks_p_values = []

    any_calls = mask_calls.sum(axis=1) > 0
    any_controls = mask_controls.sum(axis=1) > 0
    idx_calls = any_calls.index[any_calls]
    idx_controls = any_controls.index[any_controls]
    idx_inter = idx_calls.intersection(idx_controls)
    idx_union = idx_calls.union(idx_controls)

    for index in tqdm.tqdm(idx_union):
        ks_ids.append(index)
        if index in idx_inter:
            row_calls = mask_calls.loc[index]
            row_controls = mask_controls.loc[index]
            samples_pos = list(row_calls[row_calls == True].index)
            samples_neg = list(row_controls[row_controls == True].index)
            baf_pos, baf_neg = load_bafs(index, samples_pos, samples_neg, config)
        else:
            baf_pos = None
            baf_neg = None
        if baf_pos and baf_neg:
            s, p = ss.ks_2samp(sample_vals(baf_pos), sample_vals(baf_neg))
            ks_stats.append(s)
            ks_p_values.append(p)
        else:
            ks_stats.append(None)
            ks_p_values.append(None)

    result = pd.DataFrame(
        data={"baf_ks_stats": ks_stats, "baf_p_values": ks_p_values,}, index=ks_ids,
    )
    logger.info(" done performing pairwise KS tests")
    return result


def _process_deletions_par_fit(index_row, config: Config):
    index, row = index_row
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=ConvergenceWarning)
        vals = row.dropna().values
        if all(vals.shape) and vals.shape[0] >= config.gmm_n_components:
            return (
                index,
                mixture.BayesianGaussianMixture(
                    n_components=config.gmm_n_components,
                    covariance_type=config.gmm_covariance_type,
                ).fit(vals.reshape(-1, 1)),
            )
        else:
            return index, None


def process_deletions(raw_data: RawData, config: Config, caller: str = "delly"):
    logger.info("Selecting data for deletions ...")
    # ROH-positive variants.
    mask_roh = raw_data.formats["ROH"] == 1
    # Variants called by ``caller`` and not called/background without RoH
    mask_calls = raw_data.formats[caller] != 0
    mask_controls = ~mask_calls & ~mask_roh

    # DELetions.
    idx_del = raw_data.infos.SVTYPE == "DEL"
    # CNVs above reliable threshold (e.g., 5k) where read-depth signal is reliable.
    idx_cnv_reliable = raw_data.infos.SVTYPE.isin(("DEL", "DUP")) & (
        raw_data.infos.SVLEN >= config.cnv_min_reliable_size
    )
    # CNVs above a certain size are ignored (unreliable by PESR signal).
    idx_cnv_large = raw_data.infos.SVTYPE.isin(("DEL", "DUP")) & (
        raw_data.infos.SVLEN >= config.cnv_max_pesr_cnv_len
    )
    # Variants with more positive than negative samples.
    idx_too_many_called = mask_calls.sum(axis=1) > mask_calls.sum(axis=1)
    # Variants with too few controls.
    idx_too_few_controls = mask_controls.sum(axis=1) < config.min_controls

    # Select deletions of reliable size that are not too large, do not have a ROH detected, do
    # not have too many called samples, and not too few controls.
    idx = idx_del & idx_cnv_reliable & ~idx_cnv_large & ~idx_too_many_called & ~idx_too_few_controls

    # Compute log ratios for the GMM fitting.
    length = raw_data.infos.SVLEN
    thresh = np.minimum(50 / length, config.expected_snv_ratio)
    flank = np.minimum(raw_data.formats["VL"], raw_data.formats["VR"])
    middle = raw_data.formats["VM"]
    ratios = (middle.add(config.expected_snv_ratio * length, axis="rows")) / (
        flank.add(config.expected_snv_ratio * length, axis="rows")
    )
    logratios_controls = np.log10(ratios)[mask_controls][idx_del]
    ratios_calls = ratios[mask_calls][idx_del]
    logger.info("... done selecting data for deletions")

    logger.info("Fitting GMMs ...")
    with mp.Pool(config.num_procs) as p:
        gmm_models = dict(
            tqdm.tqdm(
                p.imap(
                    partial(_process_deletions_par_fit, config=config),
                    logratios_controls.iterrows(),
                ),
                total=logratios_controls.shape[0],
            ),
        )
    logger.info("... done fitting GMMs")

    logger.info("Evaluating GMMs ...")
    logratios_calls = np.log10(ratios)[mask_calls][idx_del]
    gmm_ids = []
    gmm_scores = []
    for index, row in tqdm.tqdm(logratios_calls.iterrows(), total=logratios_calls.shape[0]):
        gmm_ids.append(index)
        if gmm_models.get(index):
            vals = row.dropna().values.reshape(-1, 1)
            if all(vals.shape):
                gmm_scores.append(gmm_models[index].score(vals))
            else:
                gmm_scores.append(np.NAN)
        else:
            gmm_scores.append(np.NAN)
    logger.info("... done evaluating GMMs")

    result = pd.DataFrame(
        data={"baf_gmm_scores": gmm_scores, "baf_raw_scores": ratios_calls.mean(axis=1).values,},
        index=gmm_ids,
    )
    return result


def load_raw_data(config: Config, autosomes_only: bool = True) -> RawData:
    logger.info("Loading raw data...")

    logger.info("  Loading INFO values...")
    infos = pd.read_csv(
        config.input_pattern % "infos",
        sep="\t",
        na_values=("", ".", "nan"),
        index_col="ID",
        dtype={
            "CHROM": str,
            "POS": np.int_,
            "CHR2": str,
            "END2": np.int_,
            "ID": str,
            "SVTYPE": str,
            "SVLEN": np.int_,
        },
        nrows=config.nrows,
    )
    logger.info("  ... done loading INFO values")

    if autosomes_only:
        idx_chroms = infos.CHROM.isin(AUTOSOMES)
    else:
        idx_chroms = infos.index
    infos = infos.loc[idx_chroms]

    logger.info("  Loading FORMAT values...")
    formats = {}
    for key in config.formats:
        path = config.input_pattern % "format.{}".format(key)
        logger.info("    - %s: %s", key, path)
        formats[key] = pd.read_csv(
            path,
            index_col="ID",
            sep="\t",
            na_values=("", ".", "nan"),
            dtype={"ID": str,},
            nrows=config.nrows,
        ).loc[idx_chroms]
    logger.info("  ... done loading FORMAT values")
    logger.info("  Normalizing by autosome DoC ...")
    autosome_doc = load_autosome_doc(config)
    if set(formats[key].columns) - set(autosome_doc.keys()):
        raise Exception("no DoC info for %s" % set(formats[key].columns) - set(autosome_doc.keys()))
    for key in config.doc_normalize_formats:
        logger.info("    - %s", key)
        for col, value in autosome_doc.items():
            if col in formats[key]:
                formats[key][col] = formats[key][col] / value
    raw_data = RawData(infos=infos, formats=formats, autosome_doc=autosome_doc)
    logger.info("  ... done normalizing by autosome DoC")
    logger.info("... done loading raw data")

    if config.run_gc:
        logger.info("Running garbage collection ...")
        gc.collect()
        logger.info("  - current memory profile\n%s", hpy().heap())
        logger.info("... done running garbage collection")

    return raw_data


def create_rds_values(raw_data: RawData, config: Config, caller: str = "delly") -> pd.DataFrame:
    """Create the read depth separation values based on the median read depths."""
    logger.info("Creating RDS values ...")
    # CNVs
    idx_cnv = raw_data.infos.SVTYPE.isin(("DEL", "DUP"))
    # Size range for reliable CNV calling based on read depth.
    idx_cnv_reliable = raw_data.infos.SVLEN >= config.cnv_min_reliable_size
    # Size range that is unreliable to detect from PESR data.
    idx_cnv_large = raw_data.infos.SVLEN >= config.cnv_max_pesr_cnv_len
    # Select CNVs that can be reliable called by PESR callers and for which we can have a reliable
    # read depth signal.
    idx = idx_cnv & idx_cnv_reliable & ~idx_cnv_large

    # The calls and no-calls from the PESR ``caller``.
    mask_calls = raw_data.formats[caller] != 0

    # The read depths where ``caller`` called the variant, reliable CNVs only.
    positive_rds_median = raw_data.formats["RD"][mask_calls][idx].median(axis=1)
    negative_rds_median = raw_data.formats["RD"][~mask_calls][idx].median(axis=1)

    # Decide whether to label as positive or negative.
    rds = positive_rds_median - negative_rds_median
    result = pd.DataFrame(data={"rds_value": rds}, index=raw_data.infos.index)
    logger.info("... done creating RDS values")
    return result


def process_split_reads(raw_data, config: Config, caller: str = "delly") -> pd.DataFrame:
    def _poisson_test(row):
        """Helper function for performing Poisson test."""
        return ss.poisson.cdf(row.background, row.called)

    logger.info("Creating split read statistics ...")
    called = raw_data.formats[caller] != 0

    # For each SV, compute the median SR values of called and non-called samples.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        medians = pd.DataFrame(
            {
                "called": raw_data.formats["SR"][called].median(axis=1),
                "background": raw_data.formats["SR"][~called].median(axis=1),
            }
        )

    # Perform a Poisson test for each SV of called samples vs. background samples
    # (by median SR count) and get log-transformed P-values.
    p_values = medians.apply(_poisson_test, axis=1)
    log_pvalues = np.abs(np.log10(p_values))

    # Collect resulting feature set.
    result = pd.DataFrame(
        data={
            "sr_log_pvalue": log_pvalues,
            "sr_median_frac": medians.called / (medians.called + medians.background),
        },
        index=medians.index,
    )
    result = _fix_df(result)
    logger.info("... done creating split read statistics")
    return result


def process_rd(raw_data: RawData, config: Config, caller: str = "delly") -> pd.DataFrame:
    """... of PE/SR origin, RD origin needs limit on size..."""
    logger.info("Creating read depth statistics ...")
    idx_cnv = raw_data.infos.SVTYPE.isin(("DEL", "DUP"))
    called = raw_data.formats[caller] != 0

    rd_cnv = raw_data.formats["RD"][idx_cnv]
    rd_cnv_called = rd_cnv[called]
    rd_cnv_background = rd_cnv[~called]

    p_values = []
    # TODO: CNV variability with second-worst bin's P-value
    # second_bin_p_values = []  # TODO
    # for row_index in tqdm.tqdm(idx_cnv[idx_cnv].index):
    for row_index in idx_cnv[idx_cnv].index:
        rd_called = rd_cnv_called.loc[row_index].dropna()
        rd_background = rd_cnv_background.loc[row_index].dropna()
        # Perform power calculation to decide between two sample T-test and one-sample Z-test.
        if (
            rd_called.shape[0] <= 1
            or rd_background.shape[0] <= 1
            or (rd_background.max() == 0.0 and rd_background.min() == 0.0)
        ):
            power = -0.1
        else:
            power = sp.tt_ind_solve_power(
                nobs1=len(rd_background),
                ratio=len(rd_called) / len(rd_background),
                alpha=0.05,
                effect_size=(0.5 / rd_background.std()),
                alternative="larger",
            )
        if power >= 0.8:
            # Run two-sample T-test.
            res = ss.ttest_ind(rd_background, rd_called)
        else:
            # Run one-sample Z-test.
            rds = rd_cnv.loc[row_index].dropna()
            res = ss.ttest_1samp(rds, rds.mean())
        p_values.append(res.pvalue / 2)  # two-sided => one-sided

    # Compute log-scaled pvalues.
    log_pvalues = np.log10(p_values)

    # Compute the median separation metric.
    median_separation = rd_cnv_called.median(axis=1) - rd_cnv_background.median(axis=1)

    # Collect resulting feature set.
    result = pd.DataFrame(
        data={
            "rd_log_pvalue": log_pvalues,
            # "rd_2nd_log_pvalue": log_2nd_pvalue,  # TODO
            "rd_median_separation": median_separation,
        },
        index=rd_cnv.index,
    )
    result = _fix_df(result)
    logger.info("... done creating read depth statistics")
    return result


def process_pe(raw_data, config: Config, caller: str = "delly") -> pd.DataFrame:
    def _poisson_test(row):
        """Helper function for performing Poisson test."""
        return ss.poisson.cdf(row.background, row.called)

    logger.info("Creating PE statistics ...")
    called = raw_data.formats[caller] != 0

    # For each SV, compute the median PR values of called and non-called samples.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        medians = pd.DataFrame(
            {
                "called": raw_data.formats["PR"][called].median(axis=1),
                "background": raw_data.formats["PR"][~called].median(axis=1),
            }
        )

    # Perform a Poisson test for each SV of called samples vs. background samples
    # (by median SR count) and get log-transformed P-values.
    p_values = medians.apply(_poisson_test, axis=1)
    log_pvalues = np.abs(np.log10(p_values))

    # Collect resulting feature set.
    result = pd.DataFrame(
        data={
            "pe_log_pvalue": log_pvalues,
            "pe_median_frac": medians.called / (medians.called + medians.background),
        },
        index=medians.index,
    )
    result = _fix_df(result)
    logger.info("... done creating PE read statistics")
    return result


def process_pesr(raw_data, config: Config, caller: str = "delly") -> pd.DataFrame:
    def _poisson_test(row):
        """Helper function for performing Poisson test."""
        return ss.poisson.cdf(row.background, row.called)

    logger.info("Creating PE/SR statistics ...")
    called = raw_data.formats[caller] != 0

    # For each SV, compute the median SR values of called and non-called samples.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=RuntimeWarning)
        medians = pd.DataFrame(
            {
                "called": raw_data.formats["SR"][called].median(axis=1),
                "background": raw_data.formats["SR"][~called].median(axis=1),
            }
        ) + pd.DataFrame(
            {
                "called": raw_data.formats["PR"][called].median(axis=1),
                "background": raw_data.formats["PR"][~called].median(axis=1),
            }
        )

    # Perform a Poisson test for each SV of called samples vs. background samples
    # (by median SR count) and get log-transformed P-values.
    p_values = medians.apply(_poisson_test, axis=1)
    log_pvalues = np.abs(np.log10(p_values))

    # Collect resulting feature set.
    result = pd.DataFrame(
        data={
            "pesr_log_pvalue": log_pvalues,
            "pesr_median_frac": medians.called / (medians.called + medians.background),
        },
        index=medians.index,
    )
    result = _fix_df(result)
    logger.info("... done creating PE/SR read statistics")
    return result


def derive_baf_labels(metrics: pd.DataFrame, config: Config) -> pd.Series:
    """Derive training labels."""
    mask_fail = (config.rds_fail_min <= metrics.rds_value) & (
        metrics.rds_value < config.rds_fail_max
    )
    mask_pass = (config.rds_pass_min <= metrics.rds_value) & (
        metrics.rds_value < config.rds_pass_max
    )
    result = pd.Series(None, index=metrics.index, name="label_rds")
    result[mask_fail] = False
    result[mask_pass] = True
    return result


def derive_sr_labels(metrics: pd.DataFrame, config: Config) -> pd.Series:
    """Derive training labels."""
    mask_fail = (
        (config.rds_fail_min <= metrics.rds_value)
        & (metrics.rds_value < config.rds_fail_max)
        & (metrics.prob_baf < 0.5)
    )
    mask_pass = (
        (config.rds_pass_min <= metrics.rds_value)
        & (metrics.rds_value < config.rds_pass_max)
        & (metrics.prob_baf >= 0.9)
    )
    result = pd.Series(None, index=metrics.index, name="label_baf")
    result[mask_fail] = False
    result[mask_pass] = True
    return result


def derive_rd_labels(metrics: pd.DataFrame, config: Config) -> pd.Series:
    """Derive training labels."""
    # def label_row(self, row):
    #     if 'depth' not in row['name'] and row.svsize >= 1000:
    #         if row.BAF1_prob < 0.4 and row.SR1_prob < 0.4:
    #             return 'Fail'
    #         elif row.BAF1_prob >= 0.9 and row.SR1_prob >= 0.9:
    #             return 'Pass'
    #         else:
    #             return 'Unlabeled'
    #     elif 'depth' not in row['name'] and row.svsize < 1000:
    #         if row.SR1_prob < 0.4:
    #             return 'Fail'
    #         elif row.SR1_prob >= 0.9:
    #             return 'Pass'
    #         else:
    #             return 'Unlabeled'
    #     else:
    #         if row.BAF1_prob < 0.4:
    #             return 'Fail'
    #         elif row.BAF1_prob >= 0.9:
    #             return 'Pass'
    #         else:
    #             return 'Unlabeled'
    result = pd.Series(None, index=metrics.index, name="label_sr")
    return result


def derive_pe_labels(metrics: pd.DataFrame, config: Config) -> pd.Series:
    """Derive training labels."""
    # def label_row(self, row):
    #     if row.SR1_prob < 0.4 and row.RD_prob < 0.4 and row.svsize >= 1000:
    #         return 'Fail'
    #     elif row.SR1_prob >= 0.9 and row.RD_prob >= 0.9 and row.svsize >= 1000:
    #         return 'Pass'
    #     else:
    #         return 'Unlabeled'
    result = pd.Series(None, index=metrics.index, name="label_pe")
    return result


def derive_pesr_labels(metrics: pd.DataFrame, config: Config) -> pd.Series:
    """Derive training labels."""
    # def label_row(self, row):
    #     if (row.RD_prob < 0.4 and row.PE_prob < 0.4 and
    #             row.SR1_prob < 0.4):
    #         return 'Fail'
    #     elif (row.RD_prob >= 0.9 and row.PE_prob >= 0.9 and
    #             row.SR1_prob >= 0.9):
    #         return 'Pass'
    #     else:
    #         return 'Unlabeled'
    result = pd.Series(None, index=metrics.index, name="label_pesr")
    return result


def main():
    config = Config()
    logger.info("SV-RF")
    logger.info("config = %s", config)

    random_state = np.random.RandomState(config.random_seed)

    # Helper function for merging multiple data frames by index with OUTER JOIN.
    merge_outer = partial(pd.merge, how="outer", left_index=True, right_index=True)

    # Load the raw data.
    raw_data = load_raw_data(config)
    # Compute initial set of labels based on read-depth separation (RDS).
    rds_values = run_or_load("rds_values", config.nrows, create_rds_values, [raw_data, config])
    # Derive BAF labels from RDS values.
    baf_labels = derive_baf_labels(rds_values, config)

    # RF construction: BAF ------------------------------------------------------------------------
    logger.info("RF construction: BAF")

    # Compute BAF scores for DELetions.
    scores_del = run_or_load("baf_scores_del", config.nrows, process_deletions, [raw_data, config])
    # Compute score cutoffs for DELetion BAF scores based on initial RDS labeling.
    rf_baf_del = random_forest.SvFeatureRandomForest(
        df=scores_del.join(baf_labels),
        cols_features=["baf_gmm_scores", "baf_raw_scores"],
        col_label="label_rds",
        cols_cutoff_indeps=["baf_gmm_scores"],
        cols_cutoff_deps=["baf_raw_scores"],
    ).run()
    baf_probs_del = pd.DataFrame(data={"prob_baf": rf_baf_del.probabilites})

    # Compute BAF scores for DUPlications.
    scores_dup = run_or_load(
        "baf_scores_dup", config.nrows, process_duplications, [raw_data, config, random_state]
    )
    # Compute score cutoffs for DUPlication BAF scores based on initial RDS labeling.
    rf_baf_dup = random_forest.SvFeatureRandomForest(
        df=scores_dup.join(baf_labels),
        cols_features=["baf_p_values", "baf_ks_stats"],
        col_label="label_rds",
        cols_cutoff_indeps=["baf_p_values"],
        cols_cutoff_deps=["baf_ks_stats"],
    ).run()
    baf_probs_dup = pd.DataFrame(data={"prob_baf": rf_baf_dup.probabilites,})

    # Joint BAF-based probabilities.
    baf_probs = pd.concat([baf_probs_del, baf_probs_dup])

    # RF construction: SR -------------------------------------------------------------------------
    logger.info("RF construction: SR")

    # Perform Poisson tests on the split reads (SR) to get log-scale P-values.
    scores_sr = run_or_load("split_reads", config.nrows, process_split_reads, [raw_data, config])
    sr_labels = derive_sr_labels(reduce(merge_outer, [rds_values, baf_probs]), config)
    rf_sr = random_forest.SvFeatureRandomForest(
        df=scores_sr.join(sr_labels),
        cols_features=["sr_log_pvalue", "sr_median_frac"],
        col_label="label_baf",
        cols_cutoff_indeps=["sr_log_pvalue"],
        cols_cutoff_deps=["sr_median_frac"],
    ).run()
    logger.info("SR cutoffs are\n%s", rf_sr.cutoffs)
    sr_probs = pd.DataFrame(data={"prob_sr": rf_sr.probabilites,})

    # RF construction: RD (size-specific) ---------------------------------------------------------
    logger.info("RF construction: RD (size-specific)")

    # Load PE/SR scores
    scores_rd = run_or_load("rd_scores", config.nrows, process_rd, [raw_data, config])
    is_cnv = raw_data.infos.SVTYPE.isin(("DEL", "DUP"))
    idx_cnv = is_cnv[is_cnv].index
    above_1kbp = raw_data.infos.SVLEN >= 1000
    idx_above_1kb = above_1kbp[above_1kbp].index
    below_1kbp = raw_data.infos.SVLEN < 1000
    idx_below_1kb = below_1kbp[below_1kbp].index

    # 1. PE/SR origin CNV (>= 1kb)
    #
    # Pass: SR+BAF pass
    # Fail: SR+BAF fail
    #
    # Create labels for Pass/Fail rules.
    idx_sr_pos = sr_labels[sr_labels.label_sr == True].index
    idx_sr_neg = sr_labels[sr_labels.label_sr == False].index
    idx_baf_pos = baf_labels[baf_labels.label_baf == True].index
    idx_baf_neg = baf_labels[baf_labels.label_baf == False].index
    sr_baf_labels = pd.concat(
        [
            pd.Series(True, index=idx_sr_pos.intersection(idx_baf_pos)),
            pd.Series(False, index=idx_sr_neg.intersection(idx_baf_neg)),
        ]
    ).rename("label_sr")
    rd_labels = derive_rd_labels(config)
    # Actually do the RF based assignment
    rf_cnv_pesr_above_1kbp = random_forest.SvFeatureRandomForest(
        df=scores_rd.loc[idx_cnv & idx_above_1kb].join(rd_labels),
        cols_features=["rd_median_separation", "rd_log_pvalue"],
        col_label="label_sr",
        cols_cutoff_indeps=["rd_median_separation", "rd_log_pvalue"],
        cols_cutoff_deps=[],
    ).run()
    logger.info("RD cutoffs (PE/SR origin >= 1kb) are\n%s", rf_cnv_pesr_above_1kbp.cutoffs)
    rd_probs_pesr_above_1kb = pd.DataFrame(data={"prob_rd": rf_cnv_pesr_above_1kbp.probabilites,})

    # 2. PE/SR origin CNV (< 1kb)
    #
    # Pass: SR pass
    # Fail: SR fail
    rf_cnv_pesr_below_1kbp = random_forest.SvFeatureRandomForest(
        df=scores_rd.loc[idx_cnv & idx_below_1kb].join(rd_labels),
        cols_features=["rd_median_separation", "rd_log_pvalue"],
        col_label="label_sr",
        cols_cutoff_indeps=["rd_median_separation", "rd_log_pvalue"],
        cols_cutoff_deps=[],
    ).run()
    logger.info("RD cutoffs (PE/SR origin < 1kb) are\n%s", rf_cnv_pesr_below_1kbp.cutoffs)
    rd_probs_pesr_below_1kb = pd.DataFrame(data={"prob_rd": rf_cnv_pesr_below_1kbp.probabilites,})

    rd_probs = pd.concat([baf_probs_del, baf_probs_dup])

    # 3. Depth of Coverage origin CNV (>= 5kb)
    #
    # Pass: BAF pass
    # Fail: BAF fail
    # rf_cnv_depth_above_5kb
    # TODO: add support for CNVs from depth method
    # rf_cnv_pesr_below_1kbp = random_forest.SvFeatureRandomForest(
    #     df=scores_rd.loc[idx_cnv & ~idx_above_1kb].join(sr_labels),
    #     cols_features=["rd_median_separation", "rd_log_pvalue"],
    #     col_label="label_sr",
    #     cols_cutoff_indeps=["rd_median_separation", "rd_log_pvalue"],
    #     cols_cutoff_deps=[],
    # ).run()
    # logger.info("RD cutoffs (PE/SR origin < 1kb) are\n%s", rf_cnv_pesr_below_1kbp.cutoffs)
    # rd_labels_pesr_below_1kb = pd.DataFrame(
    #     data={
    #         "prob_rd": rf_cnv_pesr_below_1kbp.probabilites,
    #         "label_rd": (rf_cnv_pesr_below_1kbp.probabilites >= 0.5),
    #     }
    # )

    # Combine RD thresholds.
    rd_labels = pd.concat([rd_labels_pesr_above_1kb, rd_labels_pesr_below_1kb])

    # RF construction PE --------------------------------------------------------------------------
    logger.info("RF construction: PE")

    # Compute PE scores.
    score_pe = run_or_load("pe_scores", config.nrows, process_pe, [raw_data, config])
    pe_labels = derive_pe_labels(config)
    rf_pe = random_forest.SvFeatureRandomForest(
        df=score_pe.join(rd_labels),
        cols_features=["pe_log_pvalue", "pe_median_frac"],
        col_label="label_pe",
        cols_cutoff_indeps=["pe_log_pvalue"],
        cols_cutoff_deps=["pe_median_frac"],
    ).run()
    logger.info("PESR cutoffs are\n%s", rf_pe.cutoffs)
    pe_probs = pd.DataFrame(data={"prob_pe": rf_pe.probabilites})

    # RF construction: PE/SR ----------------------------------------------------------------------
    logger.info("RF construction: PE/SR")

    # Compute PE/SR scores.
    score_pesr = run_or_load("pesr_scores", config.nrows, process_pesr, [raw_data, config])
    pesr_labels = derive_pesr_labels(config)
    rf_pesr = random_forest.SvFeatureRandomForest(
        df=score_pesr.join(pesr_labels),
        cols_features=["pesr_log_pvalue", "pesr_median_frac"],
        col_label="label_rd",
        cols_cutoff_indeps=["pesr_log_pvalue"],
        cols_cutoff_deps=["pesr_median_frac"],
    ).run()
    logger.info("PESR cutoffs are\n%s", rf_pesr.cutoffs)
    pesr_probs = pd.DataFrame(data={"prob_pesr": rf_pesr.probabilites,})

    # NB: the original publication by Werling et al. (2018) describe a retraining tep of BAF and
    # SR tresholds while the code is commented out in the gatk-sv repository.

    # Consolidate Scores --------------------------------------------------------------------------
    logger.info("RF construction: score consolidation")
    merge_outer = partial(pd.merge, how="outer", left_index=True, right_index=True)
    joint_score = reduce(merge_outer, [baf_probs, sr_probs, rd_probs, pe_probs, pesr_probs])

    # Collect whether a call originated from a PESR tool.
    origin_pesr = (
        pd.DataFrame(
            {tool: raw_data.formats[tool].max(axis=1).rename(tool) for tool in config.pesr_callers}
        )
        .max(axis=1)
        .astype("bool")
    )

    # Prepare score data frame with the results.
    score = pd.DataFrame(
        {"score": None, "svtype": raw_data.infos.SVTYPE}, index=raw_data.infos.index
    )

    # Score variants of PE/SR origin <5kb.
    mask_lt5k = raw_data.infos.SVLEN < 5000
    idx_lt5k = mask_lt5k[mask_lt5k & origin_pesr].index
    score.loc[idx_lt5k, "score"] = joint_score.loc[idx_lt5k][["prob_pe", "prob_sr"]].max(axis=1)

    # Score variants of DoC origin.
    score.loc[~origin_pesr, "score"] = joint_score["prob_rd"]

    # Score variants of PE/SR origin >=5kb.
    #
    geq5kb = origin_pesr & ~mask_lt5k
    pesr_pass = (joint_score.prob_pe >= 0.5) | (joint_score.prob_sr >= 0.5)
    rd_pass = joint_score.prob_rd >= 0.5
    prob_cols = ["prob_pe", "prob_sr", "prob_rd"]

    # > If variants pass or fail both RD&PE/SR then use maximum across the three probabilities.
    score.loc[geq5kb & ~(pesr_pass ^ rd_pass), "score"] = joint_score[prob_cols].max(axis=1)

    # > If the variants pass PE/SR but not RD then use max PE/SR and reclassify as BND.
    score.loc[geq5kb & pesr_pass & ~rd_pass, "score"] = joint_score[prob_cols].max(axis=1)
    score.loc[geq5kb & pesr_pass & ~rd_pass, "svtype"] = "BND"

    # > If variants pass RD but not PE/SR then pass if they pass depth-based cutoffs.
    # >> apply rule for DELetions above 5kbp
    rescue_del = geq5kb & ~pesr_pass & rd_pass & (raw_data.infos.SVTYPE == "DEL")
    rescue_del_pass = raw_data.infos.SVLEN[rescue_del] >= 5000
    for metric, cutoff in rf_baf_del.cutoffs.items():
        rescue_del_pass = rescue_del_pass & (scores_del[metric] >= cutoff)
    score.loc[rescue_del & rescue_del_pass, "score"] = joint_score["prob_rd"]
    score.loc[rescue_del & ~rescue_del_pass, "score"] = 0.495

    # >> apply rule for DUPlications above 5kbp
    rescue_dup = geq5kb & ~pesr_pass & rd_pass & (raw_data.infos.SVTYPE == "DUP")
    rescue_dup_pass = raw_data.infos.SVLEN[rescue_dup] >= 5000
    for metric, cutoff in rf_baf_dup.cutoffs.items():
        rescue_dup_pass = rescue_dup_pass & (scores_dup[metric] >= cutoff)
    score.loc[rescue_dup & rescue_dup_pass, "score"] = joint_score["prob_rd"]
    score.loc[rescue_dup & ~rescue_dup_pass, "score"] = 0.495

    logger.info("All done. Have a nice day!")


if __name__ == "__main__":
    sys.exit(main())
