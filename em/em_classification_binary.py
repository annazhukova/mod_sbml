from collections import defaultdict
import logging
import sys
import math
import gmpy

__author__ = 'anna'


SORT_BY_PATTERN_LENGTH = 1
SORT_BY_EFM_NUMBER = 0


def ems2binary(efms, r_ids, reversible_r_ids):
    """
    Converts an EFM representation {r_id: coefficient} into a presentation {sign(coeff)*i for (r_i, coeff) \in EFM},
    that is a list of integers representing the reactions that are present in the EFM:
    for each reaction r_i participating in the EFM,
    the resulting representation will contain i if the reaction has a positive coefficient,
    -i if the reaction has a negative coefficient.
    :param efms: list of EFMs, which are represented as dictionaries {r_id: coefficient}.
    :return: list of EFMs, which are represented as sets of
    integers representing the reactions that are present in the EFM {sign(coeff)*i for (r_i, coeff) \in EFM};
             dict {r_id: i}, reaction id to the integer correspondence
    """
    int_size = math.log(sys.maxint) / math.log(2)
    result = []
    for efm in efms:
        bit_efm = []
        bit_efm_cur = 0
        i = 1
        for r_id in r_ids:
            if r_id in efm:
                coeff = efm[r_id]
            else:
                coeff = 0
            if coeff > 0:
                bit_efm_cur |= i
            bit_efm_cur, i = _increase_i(i, int_size, bit_efm, bit_efm_cur)
            if r_id in reversible_r_ids:
                if coeff < 0:
                    bit_efm_cur |= i
                bit_efm_cur, i = _increase_i(i, int_size, bit_efm, bit_efm_cur)
        bit_efm.append(bit_efm_cur)
        result.append(tuple(bit_efm))
    return result


def get_pattern_length(pattern):
    return sum(gmpy.popcount(it) for it in pattern)


def _increase_i(i, int_size, bit_efm, bit_efm_cur):
    if math.log(i) == int_size:
        bit_efm.append(bit_efm_cur)
        bit_efm_cur = 0
        i = 1
    else:
        i <<= 1
    return bit_efm_cur, i


def classify_efms(efms, r_ids, reversible_r_ids, min_pattern_len, min_efm_num=2, output_file=None, sort=SORT_BY_EFM_NUMBER):
    """
    Classifies EFMs to find common patterns.

    :param efms: list of EFMs, which are represented as a dictionaries {r_id: coefficient}.
    :param r_ids: collection of reaction ids.
    :param min_pattern_len: int, minimal length for a pattern to be considered.
    :param min_efm_num: int, at least how many EFMs should have a pattern for it to be considered.
    :param output_file: (optional) path to the file where to serialize found patterns.
    The serialization places one pattern per line. Patterns are represented as ids of the active reactions
    (for the reversed reactions, the id is preceded by minus), e.g. -R1 R3 -R7 R11 R25.
    Patterns are sorted according to the sort function.
    :param sort: an integer, representing the sorting of patterns.
    Possible values are SORT_BY_PATTERN_LENGTH and SORT_BY_EFM_NUMBER.
    :return: 2 dictionaries: p_id2efm_ids and id2pattern.
    p_id2efm_ids maps a pattern_id to ids of the EFMs containning this pattern,
    id2pattern maps a pattern_id to the pattern. Patterns are represented as sets {sign(r_i)*i for r_i \in pattern},
    that is a list of integers representing the reactions that are present in the pattern:
    for each reaction r_i participating in the pattern,
    the resulting representation will contain i if the reaction is active in its standard direction,
    -i if the reaction is reversed.
    """
    r_ids = sorted(r_ids)
    binary_ems = ems2binary(efms, r_ids, reversible_r_ids)
    logging.info("EFMs converted to %d binary vectors" % len(binary_ems))
    i = 0
    pattern2id = {}
    p_id2efm_ids = defaultdict(set)
    patterns = set()

    def process_pattern_intersection(pattern1, pattern2, new_patterns, efm_ids):
        pattern = tuple([p1 & p2 for (p1, p2) in zip(pattern1, pattern2)])
        if get_pattern_length(pattern) >= min_pattern_len:
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
    id2pattern = {p_id: p for (p, p_id) in pattern2id.iteritems()}

    if output_file:
        # serialize_pattern_hierarchy(p_id2efm_ids, i2pattern, pattern_child2parents, r_id2i, output_file)
        if SORT_BY_PATTERN_LENGTH == sort:
            sorter = lambda p_id: -get_pattern_length(id2pattern[p_id])
        elif SORT_BY_EFM_NUMBER == sort:
            sorter = lambda p_id: -len(p_id2efm_ids[p_id])
        else:
            logging.error('Unknown pattern sorting function %s. Will sort by pattern id.' % sort)
            sorter = lambda p_id: p_id
        serialize_patterns(p_id2efm_ids, id2pattern, r_ids, reversible_r_ids, output_file, sorter=sorter)

    return p_id2efm_ids, {p_id: convert_pattern(pattern, r_ids, reversible_r_ids)
                          for (p_id, pattern) in id2pattern.iteritems()}


def serialize_pattern_hierarchy(p_id2efm_ids, id2pattern, child2parents, r_id2i, output_file):
    """
    Serializes a pattern hierarchy to a file, one pattern per line. Patterns are represented as ids of the active reactions
    (for the reversed reactions, the id is preceded by minus), e.g. -R1 R3 -R7 R11 R25.
    Subpatterns are placed on the lines following their superpattern's line, preceded by a tabulation.
    :param p_id2efm_ids: dict, {pattern id: ids of EFMs containing this pattern}
    :param id2pattern: dict, {pattern id: pattern}. Patterns are presented as sets {sign(r_i)*i for r_i \in pattern},
    that is a list of integers representing the reactions that are present in the pattern:
    for each reaction r_i participating in the pattern,
    the resulting representation will contain i if the reaction is active in its standard direction,
    -i if the reaction is reversed.
    :param r_id2i: dict, {reaction_id: positive integer, associated with this reaction}
    :param output_file: path to the file where the patterns should be saved
    :param child2parents: dict, {subpattern_id: set of ids of its superpatterns}
    """
    parent2children = defaultdict(set)
    for ch_p_id, par_p_ids in child2parents.iteritems():
        for par_p_id in par_p_ids:
            parent2children[par_p_id].add(ch_p_id)
    roots = [par_p_id for par_p_id in parent2children.iterkeys() if par_p_id not in child2parents]
    i2r_id = {i: r_id for (r_id, i) in r_id2i.iteritems()}
    i2r_id.update({-i: '-%s' % r_id for (r_id, i) in r_id2i.iteritems()})

    def log_pattern(p_id, level):
        with open(output_file, 'a+') as f:
            f.write("%s%s (%d)\n" % ('\t' * level, " ".join(i2r_id[i] for i in id2pattern[p_id]),
                                     len(p_id2efm_ids[p_id])))

    dfs(roots, parent2children, level=0, process=lambda it, level: log_pattern(it, level))


def serialize_patterns(p_id2efm_ids, id2pattern, r_ids, reversible_r_ids, output_file, sorter=lambda p_id: p_id):
    """
    Serializes patterns to a file, one pattern per line. Patterns are represented as ids of the active reactions
    (for the reversed reactions, the id is preceded by minus), e.g. -R1 R3 -R7 R11 R25.
    Patterns are sorted according to the sorter function.
    :param p_id2efm_ids: dict, {pattern id: ids of EFMs containing this pattern}
    :param id2pattern: dict, {pattern id: pattern}. Patterns are presented as sets {sign(r_i)*i for r_i \in pattern},
    that is a list of integers representing the reactions that are present in the pattern:
    for each reaction r_i participating in the pattern,
    the resulting representation will contain i if the reaction is active in its standard direction,
    -i if the reaction is reversed.
    :param r_id2i: dict, {reaction_id: positive integer, associated with this reaction}
    :param output_file: path to the file where the patterns should be saved
    :param sorter: a function that given a pattern id return a value that will be used to sort a collection of patterns
    """

    def log_pattern(p_id, f):
        pattern = id2pattern[p_id]
        f.write("(%d) %s (%d)\n" % ((get_pattern_length(pattern)),
                                    ' '.join(convert_pattern(pattern, r_ids, reversible_r_ids)),
                                    len(p_id2efm_ids[p_id])))

    with open(output_file, 'w+') as f:
        for p_id in sorted(p_id2efm_ids.iterkeys(), key=sorter):
            log_pattern(p_id, f)


def convert_pattern(pattern, r_ids, reversible_r_ids):
    converted_pattern = []
    i = 1
    pattern_iterable = iter(pattern)
    cur_pattern = next(pattern_iterable)
    int_size = math.log(sys.maxint) / math.log(2)
    for r_id in r_ids:
        if cur_pattern & i:
            converted_pattern.append(r_id)
        if math.log(i) == int_size:
            cur_pattern = next(pattern_iterable)
            i = 1
        else:
            i <<= 1
        if r_id in reversible_r_ids:
            if cur_pattern & i:
                converted_pattern.append('-' + r_id)
            if math.log(i) == int_size:
                cur_pattern = next(pattern_iterable)
                i = 1
            else:
                i <<= 1
    return converted_pattern


def dfs(roots, element2children, process, level):
    """
    Preorder depth-first search.
    :param roots: an iterable of root elements
    :param element2children: dict, {element: iterable of child elements}
    :param process: the operation to be performed on each element, lambda element, level: operation(element, level)
    :param level: the level of the current root in the complete tree.
    """
    for it in roots:
        process(it, level)
        if it in element2children:
            dfs(element2children[it], element2children, process, level + 1)


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




