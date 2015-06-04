from collections import defaultdict

from em.efm_manager import binary2efm, get_int_size

__author__ = 'anna'


def classify_reactions_by_efm(id2efm, r_ids, rev_r_ids):
    """
    Calculates the reaction id to EFM ids correspondence.

    :param id2efm: dictionary that maps EFM identifiers (int) to the EFMs,
    which are represented as tuples (binary_representation, non_zero_coefficients).

    A binary representation of an EFM is a list of integers whose binary representations
    correspond to the reactions that are active in the EFM: if the reaction is active,
    the corresponding bit is set to 1.
    If the total number of reactions in the model is larger that the number of bits in an int, several ints are used.

    Example: For a model containing reactions r1, r2(reversible), r3(reversible), r4, r5,
    a EFM: 3 r1, -2 r2, 1 r3, 1 r5 would be represented as [77], as the binary representation of 77 is '1001101'
    that corresponds to '1 r5, 0 r4, 0 -r3, 1 r3, 1 -r2, 0 r2, 1 r1', where -ri is the reversed version of ri.
    The non_zero_coefficients for this EFM are [3, -2, 1, 1].

    :return: dict: r_id: efm_ids, where the reversed reactions are represented as "-%s" % r_id.
    """
    int_size = get_int_size()
    id2efm = {efm_id: efm[0] for (efm_id, efm) in id2efm.iteritems()}
    r_id2efm_ids = defaultdict(set)
    for efm_id, efm in id2efm.iteritems():
        r_id2coeff = binary2efm(efm, r_ids, rev_r_ids, int_size=int_size, binary=True)
        for r_id, coeff in r_id2coeff.iteritems():
            if r_id in rev_r_ids and coeff < 0:
                r_id2efm_ids['-' + r_id].add(efm_id)
            else:
                r_id2efm_ids[r_id].add(efm_id)
    return r_id2efm_ids