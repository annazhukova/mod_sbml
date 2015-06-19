from collections import defaultdict
import logging
from constraint_based_analysis.efm.EFM import EFM, get_binary_efm_len

__author__ = 'anna'


def classify_efms(id2efm, min_pattern_len, min_efm_num=2):
    """
    Classifies EFMs to find common patterns.

    :param efms: list of EFMs, which are represented as tuples (binary_representation, non_zero_coefficients).

    A binary representation of an EFM is a list of integers whose binary representations
    correspond to the reactions that are active in the EFM: if the reaction is active,
    the corresponding bit is set to 1.
    If the total number of reactions in the model is larger that the number of bits in an int, several ints are used.

    Example: For a model containing reactions r1, r2(reversible), r3(reversible), r4, r5,
    a EFM: 3 r1, -2 r2, 1 r3, 1 r5 would be represented as [77], as the binary representation of 77 is '1001101'
    that corresponds to '1 r5, 0 r4, 0 -r3, 1 r3, 1 -r2, 0 r2, 1 r1', where -ri is the reversed version of ri.
    The non_zero_coefficients for this EFM are [3, -2, 1, 1].

    :param min_pattern_len: int, minimal length for a pattern to be considered.

    :param min_efm_num: int, at least how many EFMs should have a pattern for it to be considered.

    :return: 2 dictionaries: p_id2efm_ids and id2pattern.
    p_id2efm_ids maps a pattern_id to ids of the EFMs containing this pattern,
    id2pattern maps a pattern_id to the binary representation of the pattern.

    A binary representation of a pattern is a list of integers whose binary representations
    correspond to the reactions that are active in this pattern: if the reaction is active,
    the corresponding bit is set to 1.
    If the total number of reactions in the model is larger that the number of bits in an int, several ints are used.

    Example: For a model containing reactions r1, r2(reversible), r3(reversible), r4, r5,
    a pattern: r1, r2, r3, r5 would be represented as [77], as the binary representation of 77 is '1001101'
    that corresponds to '1 r5, 0 r4, 0 -r3, 1 r3, 1 -r2, 0 r2, 1 r1', where -ri is the reversed version of ri.
    """
    id2binary_ems = {efm_id: efm.binary_efm for (efm_id, efm) in id2efm.iteritems()}
    binary_ids = sorted(id2binary_ems.iterkeys())
    logging.info("EFMs converted to %d binary vectors" % len(id2binary_ems))
    pattern2id = {}
    p_id2efm_ids = defaultdict(set)
    patterns = set()

    def process_pattern_intersection(pattern1, pattern2, new_patterns, efm_ids):
        pattern = tuple([p1 & p2 for (p1, p2) in zip(pattern1, pattern2)])
        if get_binary_efm_len(pattern) >= min_pattern_len:
            if pattern not in pattern2id:
                new_patterns.add(pattern)
                pattern2id[pattern] = len(pattern2id)
            p_id = pattern2id[pattern]
            p_id2efm_ids[p_id] |= efm_ids
            return p_id
        return None

    logging.info("Calculating patterns contained in at least 2 EFMs...")
    i = 0
    for efm_id1 in binary_ids:
        efm1 = id2binary_ems[efm_id1]
        for efm_id2 in binary_ids[i + 1:]:
            efm2 = id2binary_ems[efm_id2]
            process_pattern_intersection(efm1, efm2, patterns, {efm_id1, efm_id2})
        i += 1
    logging.info("... found %d patterns." % len(patterns))
    level = 3
    while patterns:
        logging.info("Calculating subpatterns contained in at least %d EFMs..." % level)
        new_patterns = set()
        for pattern in patterns:
            p_id = pattern2id[pattern]
            pattern = set(pattern)
            efm_ids = p_id2efm_ids[p_id]
            for efm_id in binary_ids:
                if efm_id not in efm_ids:
                    process_pattern_intersection(id2binary_ems[efm_id], pattern, new_patterns, {efm_id} | efm_ids)
        patterns = new_patterns
        logging.info("... found %d patterns." % len(patterns))
        level += 1
    id2pattern = {p_id: p for (p, p_id) in pattern2id.iteritems() if len(p_id2efm_ids[p_id]) >= min_efm_num}
    id2id = dict(zip(id2pattern.iterkeys(), xrange(1, len(id2pattern) + 1)))
    sample_efm = next(id2efm.itervalues())
    id2pattern = {id2id[p_id]: EFM(r_ids=sample_efm.r_ids, rev_r_ids=sample_efm.rev_r_ids, int_size=sample_efm.int_size,
                                   binary_efm=p) for (p_id, p) in id2pattern.iteritems()}
    p_id2efm_ids = {id2id[p_id]: efm_ids for (p_id, efm_ids) in p_id2efm_ids.iteritems() if p_id in id2id}
    return p_id2efm_ids, id2pattern


def filter_patterns(min_efm_num, p_id2efm_ids, pattern2id):
    """
    Filters patterns to eliminate those that are contained in less than min_efm_num EFMS.
    :param min_efm_num: minimal number of EFMs that need to contain a pattern for it to be kept
    :param p_id2efm_ids: dict, {pattern_id : set of ids of EFMs that contain pattern with this pattern_id}
    :param pattern2id: dict, {pattern: pattern_id}
    """
    logging.info("Filtering found patterns to eliminate those that are contained in less than %d EFMs..." % min_efm_num)
    patterns = set(pattern2id.keys())
    for pattern in patterns:
        p_id = pattern2id[pattern]
        if len(p_id2efm_ids[p_id]) < min_efm_num:
            del p_id2efm_ids[p_id]
            del pattern2id[pattern]




