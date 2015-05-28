import logging
import os
import math
import sys
from em.em_manager import binary2efm

__author__ = 'anna'


def efms2acom_input(efms, r_ids, rev_r_ids, react_file, binary_efm_file):
    with open(react_file, 'w+') as f:
        f.write(','.join(r_ids) + '\n')
    logging.info('Wrote reactions to %s.' % react_file)
    int_size = math.log(sys.maxint) / math.log(2)
    with open(binary_efm_file, 'w+') as f:
        for (binary_efm, coefficients) in efms:
            r_id2coefficient = binary2efm(binary_efm, r_ids, rev_r_ids, int_size, coefficients, binary=True)
            f.write(' '.join((str(r_id2coefficient[r_id] if r_id in r_id2coefficient else 0) for r_id in r_ids)) + '\n')
    logging.info('Wrote EFMs matrix to %s.' % binary_efm_file)


def acom_classification(efms, r_ids, rev_r_ids, directory, similarity_threshold, min_pattern_size,
                        acom_path='/home/anna/Applications/acom-c/acom-c'):
    react_file = '%s/acom_reactions.react' % directory
    efm_file = '%s/acom_efms.mat' % directory
    efms2acom_input(efms, r_ids, rev_r_ids, react_file, efm_file)
    os.system("%s %s %s %d %d %d" %
              (acom_path, react_file, efm_file, len(efms), similarity_threshold, min_pattern_size))
