import logging

from mod_sbml.constraint_based_analysis.efm.efm_manager import compute_efms

TREEEFM_PATH = "/home/anna/Applications/TreeEFM/tool/TreeEFMseq"

__author__ = 'anna'

ZERO_THRESHOLD = 1e-6

def get_efms(target_r_id, target_r_reversed, r_id2rev, sbml, directory, max_efm_number=1000, threshold=ZERO_THRESHOLD,
             efms=None, r_ids=None, rewrite=True, tree_efm_path=TREEEFM_PATH):
    if not efms:
        if tree_efm_path:
            efms, r_ids, rev_r_ids = compute_efms(sbml, directory, max_efm_number, target_r_id, target_r_reversed,
                                                  tree_efm_path, r_id2rev, threshold=threshold, r_ids=r_ids,
                                                  rewrite=rewrite)
            if not efms:
                logging.info('Found no EFMs of interest.')
                return None
        else:
            raise ValueError('No EFMs found :(. Probably, you forgot to specify TreeEFM path?')
    id2efm = dict(zip(xrange(1, len(efms) + 1), efms))
    return id2efm


def calculate_imp_rn_threshold(all_id2efm, essential_rn_number, imp_rn_threshold):
    if imp_rn_threshold is None:
        imp_rn_threshold = int(len(all_id2efm) * 0.3)
        if essential_rn_number:
            at_least_num = min(int(essential_rn_number * 1.25), len(all_id2efm))
            if at_least_num == len(all_id2efm):
                at_least_num = min(max(len(all_id2efm) / 2, essential_rn_number + 2), len(all_id2efm) - 1)
            imp_rn_threshold = max(2, at_least_num, imp_rn_threshold)
    return imp_rn_threshold


def calculate_min_pattern_len(avg_efm_len, essential_rn_number, min_efm_len, min_pattern_len):
    if not min_pattern_len:
        min_pattern_len = min(int(avg_efm_len * 0.3), min_efm_len)
        if essential_rn_number:
            min_pattern_len = min(max(essential_rn_number + 3, min_pattern_len), avg_efm_len)
    return min_pattern_len


def calculate_min_clique_len(avg_efm_len, essential_rn_number, min_clique_len, min_efm_len):
    if not min_clique_len:
        min_clique_len = min(avg_efm_len / 4, min_efm_len)
        if essential_rn_number:
            min_clique_len = min(max(essential_rn_number + 3, min_clique_len), avg_efm_len)
    return min_clique_len


