import numpy as np


def result_reader(samples, scores, m_factor=20, nratio=1., seed=None):
    """
    Calculate the index of accepted proposal draws based on the scores and jacobians calcualted above
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

    return inds_field, inds_clust, inds_2d


def result_reader2(samples, scores, m_factor=20, nratio=1., seed=None):

    dc_score = np.exp(scores["dc"]) * np.abs(scores["dc_jac"])
    wr_score = np.exp(scores["wr"]) * np.abs(scores["wr_jac"])
    wcr_clust_score = np.exp(scores["wcr_clust"]) * np.abs(scores["wcr_clust_jac"])
    wcr_rands_score = np.exp(scores["wcr_rands"]) * np.abs(scores["wcr_rands_jac"])
    #     fcls = fclfunc(samples["LOGR"]) - 1.

    #     umult = np.max((fcls, np.ones(len(fcls))), axis=0)
    #     print(umult)
    rng = np.random.RandomState(seed)
    uniform = rng.uniform(0, 1, len(samples))  # * umult

    p_proposal = m_factor * dc_score * wr_score
    p_rands_ref = wcr_rands_score
    p_clust_ref = wcr_clust_score  # * fcls
    print(nratio)
    inds_field = (uniform < (p_rands_ref / nratio / p_proposal))
    inds_clust = ((p_rands_ref / nratio / p_proposal) < uniform) * (uniform < (p_clust_ref  / p_proposal))
    inds_2d = ((uniform < (p_clust_ref  / p_proposal)))

    return inds_field, inds_clust, inds_2d

