from collections import Counter
import logging

from networkx import find_cliques, Graph
from mod_sbml.constraint_based_analysis.efm.EFM import EFM

__author__ = 'anna'


def detect_patterns(id2efm, min_pattern_len, efm_num=2):
    """
    Calculates the reaction id to EFM ids correspondence.

    :param id2efm: dictionary that maps EFM identifiers (int) to the EFMs.

    :return: dict: r_id: efm_ids, where the reversed reactions are represented as "-%s" % r_id.
    """
    logging.info("Going to rank reactions by EFM number.")
    r_id_pair2count = Counter()

    for efm_id, efm in id2efm.iteritems():
        r_id2coeff = efm.to_r_id2coeff(binary=True)
        r_ids = [r_id if coeff > 0 else '-' + r_id for (r_id, coeff) in r_id2coeff.iteritems()]
        i = 0
        for r_id in r_ids:
            i += 1
            for r_id2 in r_ids[i:]:
                r_id_pair2count.update({tuple(sorted([r_id, r_id2])): 1})
    gr = Graph()
    for (r_id1, r_id2), count in r_id_pair2count.iteritems():
        if count >= efm_num:
            gr.add_edge(r_id1, r_id2)
    sample_efm = next(id2efm.itervalues())
    r_ids, rev_r_ids, int_size = sample_efm.r_ids, sample_efm.rev_r_ids, sample_efm.int_size
    clique2r_id2coeff = lambda clique: {(r_id if '-' != r_id[0] else r_id[1:]): (1 if '-' != r_id[0] else -1)
                                        for r_id in clique}
    patterns = [EFM(r_ids=r_ids, rev_r_ids=rev_r_ids, int_size=int_size,
                    r_id2coeff=clique2r_id2coeff(clique)) for clique in
                (clique for clique in find_cliques(gr) if len(clique) >= min_pattern_len)]
    id2pattern = dict(zip(xrange(1, len(patterns) + 1), patterns))
    p_id2efm_ids = {p_id: {efm_id for efm_id in id2efm.iterkeys()
                           if p.intersection(id2efm[efm_id]) == p} for (p_id, p) in id2pattern.iteritems()}
    all_efm_intersection = reduce(lambda p1, p2: p1.intersection(p2), id2efm.itervalues(), next(id2efm.itervalues()))
    return p_id2efm_ids, id2pattern, all_efm_intersection




