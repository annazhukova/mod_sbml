import logging
import os
import math
import sys
from em.efm_manager import binary2efm

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
    """
    Classifies binary elementary flux modes (EFMs) using ACoM method:
    Peres et al. 2011 (doi:10.1016/j.biosystems.2010.12.001).
    :param efms: a list of EFMs in binary form.

    A binary representation of an EFM is a list of integers whose binary representations
    correspond to the reactions that are active in the EFM: if the reaction is active,
    the corresponding bit is set to 1.
    If the total number of reactions in the model is larger that the number of bits in an int, several ints are used.

    Example: For a model containing reactions r1, r2(reversible), r3(reversible), r4, r5,
    a EFM: 3 r1, -2 r2, 1 r3, 1 r5 would be represented as [77], as the binary representation of 77 is '1001101'
    that corresponds to '1 r5, 0 r4, 0 -r3, 1 r3, 1 -r2, 0 r2, 1 r1', where -ri is the reversed version of ri.

    :param r_ids: ordered collection of reaction ids (strings).

    :param reversible_r_ids: set of ids of reversible reactions (strings).

    :param directory: path to the directory where the intermediate files needed foe ACoM will be saved.

    :param similarity_threshold: int, at least how many common elements two EFMs should have
    to be considered as neighbours (see Peres et al. 2011).

    :param min_pattern_size: int, minimal motif length to be considered.

    :param acom_path: path to the acom-c programme
    """
    react_file = '%s/acom_reactions.react' % directory
    efm_file = '%s/acom_efms.mat' % directory
    efms2acom_input(efms, r_ids, rev_r_ids, react_file, efm_file)
    os.system("%s %s %s %d %d %d" %
              (acom_path, react_file, efm_file, len(efms), similarity_threshold, min_pattern_size))
