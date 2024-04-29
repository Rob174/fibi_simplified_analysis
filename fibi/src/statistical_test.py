"""
Contains all logic related to statistical tests: 
- mapping the effect size and pvalues to their respective categories, 
- run the wiloconx test
- is_conclusion_sign_verified depending of if we maximize or minimize the objective
"""
import numpy as np
from typing import *
from scipy.stats import distributions
from scipy.stats._stats_py import rankdata, find_repeats

if TYPE_CHECKING:
    effect_size_category = Literal["nan","small","small-medium","medium-big","big"]
    pvalue_category = Literal["nan","small","big"]
    test_category = Literal["ZTest", "TTest", "Wilcoxon", "SignTest"]

def mapping_effect_size_category(value: float, test: 'test_category') -> 'effect_size_category':
    """Given the effect size value and the test, gives if the effect size is small, small-medium ..."""
    if np.isnan(value):
        return "nan"
    value = abs(value)
    if test in ["ZTest", "TTest"]:
        if value < 0.2:
            return "small"
        if 0.2 <= value < 0.5:
            return "small-medium"
        elif 0.5 <= value < 0.8:
            return "medium-big"
        else:
            return "big"
    elif test in ["Wilcoxon", "SignTest"]:
        if value < 0.1:
            return "small"
        if 0.1 <= value < 0.3:
            return "small-medium"
        elif 0.3 <= value < 0.5:
            return "medium-big"
        else:
            return "big"
    else:
        raise Exception("Unknown test " + test)

def mapping_pvalue_category(value: float) -> 'pvalue_category':
    """Given the pvalue, gives if the pvalue is small, big or nan"""
    if np.isnan(value):
        return "nan"
    if value < 0.05:
        return "small"
    else:
        return "big"

def is_conclusion_sign_verified(avg: float, maximization: bool = True, init_random: bool = True) -> bool:
    """Given the average and wether if we try to maximize the objective, tells if the conclusion of Hansen is verified
    If initialization is random then FI gives better results
    -> for maximization means avg((BI-FI)/init) < 0
    -> for minimization means avg((BI-FI)/init) > 0
    If initialization is greedy then BI gives better results
    -> for maximization means avg((BI-FI)/init) > 0
    -> for minimization means avg((BI-FI)/init) < 0
    """
    # (maximization and init_random) XOR (not maximization and not init_random)
    if maximization == init_random:
        mapping = {
            -1: True,
            1: False,
            0: False
        }
    else:
        mapping = {
            1: True,
            -1: False, 
            0: False 
        }
    return mapping[int(np.sign(avg))]

def run_wilcoxon(diff: np.ndarray) -> Dict[str, Any]:
    """Run the wilcoxon test and returns a dictionnary with especially the pvalue and effect_size"""
    d = diff.astype(np.float64)  # type: ignore
    assert len(diff) >= 10, f"Warning: sample size too small {len(diff)}"

    mode = "approx"

    count = len(d)
    n_zero = np.sum(d == 0)
    if n_zero == len(d):
        prob = np.nan
        return dict(
            pvalue=np.nan,
            effect_size=np.nan,
            count=count,
            n_zero=len(d),
            r_plus=0,
            r_minus=0,
            T=0,
            mode=mode,
            prob=np.nan,
            mn=None,
            se=None,
            z=None,
        )

    r = rankdata(abs(d))
    r_plus = np.sum((d > 0) * r)
    r_minus = np.sum((d < 0) * r)

    # return min for two-sided test, but r_plus for one-sided test
    # the literature is not consistent here
    # r_plus is more informative since r_plus + r_minus = count*(count+1)/2,
    # i.e. the sum of the ranks, so r_minus and the min can be inferred
    # (If alternative='pratt', r_plus + r_minus = count*(count+1)/2 - r_zero.)
    # [3] uses the r_plus for the one-sided test, keep min for two-sided test
    # to keep backwards compatibility

    T = min(r_plus, r_minus)
    mn = None
    se = None
    z = None
    rep = None
    mn = count * (count + 1.0) * 0.25
    se = count * (count + 1.0) * (2.0 * count + 1.0)

    replist, repnum = find_repeats(r)
    if repnum.size != 0:
        rep = repnum
        # Correction for repeated elements.
        se -= 0.5 * (repnum * (repnum * repnum - 1)).sum()

    se = np.sqrt(se / 24)

    # apply continuity correction if applicable
    d = 0

    # compute statistic and p-value using normal approximation
    z = (T - mn - d) / se
    prob = 2.0 * distributions.norm.sf(abs(z))
    es = 4 * abs(T - (r_plus - r_minus) / 2) / (count * (count + 1))
    return dict(
        effect_size=es,
        pvalue=prob,
        n_zero=n_zero,
        r_plus=r_plus,
        r_minus=r_minus,
        T=T,
        mode=mode,
        prob=prob,
        mn=mn,
        se=se,
        rep=rep,
        z=z
    )
    