from collections import Counter
import logging

from networkx import find_cliques, Graph
from mod_sbml.constraint_based_analysis.efm.EFM import EFM

__author__ = 'anna'


def detect_cliques(id2efm, min_clique_size, efm_num=2):
    """
    The method takes the found EFMs and constructs the reaction graph in the following way:
    nodes are marked with reaction ids (or -r_id for the reversed versions of reversible reactions),
    there exists an edge between nodes r_i and r_j iff they are related, i.e.
    there exists at least efm_num of EFMs that contain both reaction r_i and reaction r_j.
    The method then detects the maximal cliques of size greater or equal to min_clique_size.

    :param id2efm: dictionary that maps EFM identifiers (int) to the EFMs.
    :param min_clique_size: int, minimal size of a clique for it to be considered.
    :param efm_num: int (optional, default value is 2), minimal number of EFMs that should contain two reactions
    for them to be considered related.

    :return: id2clique
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
    cliques = [EFM(r_ids=r_ids, rev_r_ids=rev_r_ids, int_size=int_size,
                    r_id2coeff=clique2r_id2coeff(clique)) for clique in
                (clique for clique in find_cliques(gr) if len(clique) >= min_clique_size)]
    id2clique = dict(zip(xrange(1, len(cliques) + 1), cliques))
    return id2clique




