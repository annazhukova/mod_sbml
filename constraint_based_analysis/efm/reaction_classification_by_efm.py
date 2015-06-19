from collections import defaultdict
import logging

__author__ = 'anna'


def classify_reactions_by_efm(id2efm):
    """
    Calculates the reaction id to EFM ids correspondence.

    :param id2efm: dictionary that maps EFM identifiers (int) to the EFMs.

    :return: dict: r_id: efm_ids, where the reversed reactions are represented as "-%s" % r_id.
    """
    logging.info("Going to rank reactions by EFM number.")
    r_id2efm_ids = defaultdict(set)
    for efm_id, efm in id2efm.iteritems():
        r_id2coeff = efm.to_r_id2coeff(binary=True)
        for r_id, coeff in r_id2coeff.iteritems():
            if coeff < 0:
                r_id2efm_ids['-' + r_id].add(efm_id)
            else:
                r_id2efm_ids[r_id].add(efm_id)
    return r_id2efm_ids