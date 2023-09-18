import numpy as np


def result_reader(samples, scores, m_factor=20, nratio=1., seed=None):
    """
    Calculate the index of accepted proposal draws based on the scores and jacobians

    Parameters
    ----------
    samples: pd.DataFrame
        random samples from the rejection sampling
    scores: pd.DataFrame
        log probabilities and related quantities from the rejection sampling
    m_factor: float
        We are going to accept 1 / m_factor fraction  of the points, that is with 20, we accept 5%
    nratio: float
        overdensity of points around clusters compared to random points
    seed: int
        random seed

    Returns
    -------
        Indexes of samples assigned to the field, Indexes of samples assigned to clusters


    """
    dc_score = np.exp(scores["dc"]) * np.abs(scores["dc_jac"])
    wr_score = np.exp(scores["wr"]) * np.abs(scores["wr_jac"])
    wcr_clust_score = np.exp(scores["wcr_clust"]) * np.abs(scores["wcr_clust_jac"])
    wcr_rands_score = np.exp(scores["wcr_rands"]) * np.abs(scores["wcr_rands_jac"])

    rng = np.random.RandomState(seed)
    uniform = rng.uniform(0, 1, len(samples))

    p_proposal = m_factor * dc_score * wr_score
    p_rands_ref = wcr_rands_score
    p_clust_ref = wcr_clust_score
    print(nratio)
    inds_field = (uniform < (p_rands_ref / nratio / p_proposal))
    inds_clust = ((p_rands_ref / nratio / p_proposal) < uniform) * (uniform < (p_clust_ref  / p_proposal))
    inds_2d = ((uniform < (p_clust_ref  / p_proposal)))

    return inds_field, inds_clust,