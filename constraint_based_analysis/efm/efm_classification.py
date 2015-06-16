from collections import defaultdict
import logging

from constraint_based_analysis.efm.efm_manager import get_binary_efm_length

__author__ = 'anna'


def classify_efms(efms, min_pattern_len, min_efm_num=2):
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
    binary_ems = [binary_efm for (binary_efm, coeffs) in efms]
    logging.info("EFMs converted to %d binary vectors" % len(binary_ems))
    i = 0
    pattern2id = {}
    p_id2efm_ids = defaultdict(set)
    patterns = set()

    def process_pattern_intersection(pattern1, pattern2, new_patterns, efm_ids):
        pattern = tuple([p1 & p2 for (p1, p2) in zip(pattern1, pattern2)])
        if get_binary_efm_length(pattern) >= min_pattern_len:
            if pattern not in pattern2id:
                new_patterns.add(pattern)
                pattern2id[pattern] = len(pattern2id)
            p_id = pattern2id[pattern]
            p_id2efm_ids[p_id] |= efm_ids
            return p_id
        return None

    logging.info("Calculating patterns contained in at least 2 EFMs...")
    for efm1 in binary_ems:
        j = i + 1
        for efm2 in binary_ems[i + 1:]:
            process_pattern_intersection(efm1, efm2, patterns, {i, j})
            j += 1
        i += 1
    logging.info("... found %d patterns." % len(patterns))
    level = 3
    pattern_child2parents = defaultdict(set)
    while patterns:
        logging.info("Calculating subpatterns contained in at least %d EFMs..." % level)
        new_patterns = set()
        for pattern in patterns:
            p_id = pattern2id[pattern]
            pattern = set(pattern)
            efm_ids = p_id2efm_ids[p_id]
            i = 0
            for efm in binary_ems:
                if i not in efm_ids:
                    pp_id = process_pattern_intersection(efm, pattern, new_patterns, {i} | efm_ids)
                    if pp_id and pp_id != p_id:
                        pattern_child2parents[pp_id].add(p_id)
                i += 1
        patterns = new_patterns
        logging.info("... found %d patterns." % len(patterns))
        level += 1
    filter_patterns(min_efm_num, p_id2efm_ids, pattern2id, pattern_child2parents)
    id2id = dict(zip(pattern2id.itervalues(), xrange(1, len(pattern2id) + 1)))
    id2pattern = {id2id[p_id]: p for (p, p_id) in pattern2id.iteritems()}
    p_id2efm_ids = {id2id[p_id]: efm_ids for (p_id, efm_ids) in p_id2efm_ids.iteritems()}
    return p_id2efm_ids, id2pattern


def filter_patterns(min_efm_num, p_id2efm_ids, pattern2id, pattern_child2parents):
    """
    Filters patterns to eliminate those that are contained in less than min_efm_num EFMS.
    :param min_efm_num: minimal number of EFMs that need to contain a pattern for it to be kept
    :param p_id2efm_ids: dict, {pattern_id : set of ids of EFMs that contain pattern with this pattern_id}
    :param pattern2id: dict, {pattern: pattern_id}
    :param pattern_child2parents: dict, {subpattern_id: set of ids of its superpatterns}
    """
    logging.info("Filtering found patterns to eliminate those that are contained in less than %d EFMs..." % min_efm_num)
    patterns = set(pattern2id.keys())
    for pattern in patterns:
        p_id = pattern2id[pattern]
        if len(p_id2efm_ids[p_id]) < min_efm_num:
            del p_id2efm_ids[p_id]
            del pattern2id[pattern]
            if p_id in pattern_child2parents:
                del pattern_child2parents[p_id]
    filtered_p_ids = set(pattern2id.itervalues())
    for ch_p_id in filtered_p_ids:
        if ch_p_id in pattern_child2parents:
            pattern_child2parents[ch_p_id] &= filtered_p_ids
            if not pattern_child2parents[ch_p_id]:
                del pattern_child2parents[ch_p_id]




