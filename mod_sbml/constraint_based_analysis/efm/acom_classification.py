from dircache import listdir
import logging
import os
import re

from mod_sbml.constraint_based_analysis.efm import efm_manager

__author__ = 'anna'


def efms2acom_input(efm_id2efms, r_ids, react_file, binary_efm_file):
    with open(react_file, 'w+') as f:
        f.write(','.join(r_ids) + '\n')
    logging.info('Wrote reactions to %s.' % react_file)

    with open(binary_efm_file, 'w+') as f:
        for efm_id, efm in sorted(efm_id2efms.iteritems(), key=lambda (efm_id, _): efm_id):
            r_id2coefficient = efm.to_r_id2coeff(binary=True)
            f.write(' '.join((str(r_id2coefficient[r_id] if r_id in r_id2coefficient else 0) for r_id in r_ids)) + '\n')
    logging.info('Wrote EFMs matrix to %s.' % binary_efm_file)



def acom_classification(efm_id2efms, r_ids, directory, similarity_threshold, min_pattern_length,
                        acom_path='/home/anna/Applications/acom-c/acom-c'):
    """
    Classifies binary elementary flux modes (EFMs) using ACoM method:
    Peres et al. 2011 (doi:10.1016/j.biosystems.2010.12.001).
    :param efm_id2efms: dictionary that maps efm_id (int) to EFM.

    :param r_ids: ordered collection of reaction ids (strings).

    :param directory: path to the directory where the intermediate files needed foe ACoM will be saved.

    :param similarity_threshold: int, at least how many common elements two EFMs should have
    to be considered as neighbours (see Peres et al. 2011).

    :param min_pattern_length: int, minimal motif length to be considered.

    :param acom_path: path to the acom-c programme
    """
    logging.info("Going to classify EFMs with ACoM: min motif length is set to %d, similarity threshold is set to %d."
                 % (min_pattern_length, similarity_threshold))
    react_file = os.path.join(directory, 'acom_reactions.react')
    efm_file = os.path.join(directory, 'acom_efms.mat')
    efms2acom_input(efm_id2efms, r_ids, react_file, efm_file)
    os.system("%s %s %s %d %d %d" %
              (acom_path, react_file, efm_file, len(efm_id2efms), similarity_threshold, min_pattern_length))
    process_clusters(directory, efm_id2efms, similarity_threshold, min_pattern_length)


def process_clusters(dest_path, efm_id2efm, similarity_threshold, min_pattern_size):
    classes_dir = os.path.dirname(os.path.abspath(efm_manager.__file__))
    for f in sorted(listdir(classes_dir)):
        if re.search('^class\d+.txt$', f):
            os.rename(os.path.join(classes_dir, f), os.path.join(dest_path, f))
    logging.info('Moved ACoM classes files to %s.' % dest_path)

    id2pattern = {}
    p_id2efm_ids = {}
    found_efm_ids = set()
    for f in sorted(listdir(dest_path)):
        if re.search('^class\d+.txt$', f):
            pattern_id = int(next(re.finditer('\d+', f)).group(0))
            efm_ids = set()
            pattern = None
            with open(os.path.join(dest_path, f), 'r') as f:
                for line in f:
                    if line.startswith('Class'):
                        continue
                    line = line.replace('\n', '').strip()
                    if not line:
                        continue
                    efm_id, efm = line.split(': ')
                    efm_id = int(efm_id)
                    efm_r_ids = efm.split(' ')
                    efm_ids.add(efm_id)
                    if not pattern:
                        pattern = set(efm_r_ids)
                    else:
                        pattern &= set(efm_r_ids)
            id2pattern[pattern_id] = pattern
            p_id2efm_ids[pattern_id] = efm_ids
            found_efm_ids |= efm_ids
    logging.info('Processed ACoM classes.')

    outlier_efm_ids = set(efm_id2efm.iterkeys()) - found_efm_ids
    if outlier_efm_ids and len(outlier_efm_ids) < len(efm_id2efm):
        outlier_file = os.path.join(dest_path, 'outliers.txt')
        with open(outlier_file, 'w+') as f:
            f.write('Outliers:\n')
            for efm_id in outlier_efm_ids:
                f.write('%d: %s\n' % (efm_id, efm_id2efm[efm_id].to_string(binary=True)))
        logging.info('Saved outliers to %s.' % outlier_file)

    motif_file = os.path.join(dest_path, 'motifs.txt')
    with open(motif_file, 'w+') as f:
        f.write('Performed ACoM using %d as similarity threshold and %d as min motif length.\n\n'
                % (similarity_threshold, min_pattern_size))
        if id2pattern:
            f.write('Found %d motifs:\n' % len(id2pattern))
            for p_id, pattern in id2pattern.iteritems():
                f.write('%d\t(len=%d,\tnum of EFMs=%d):\t%s\n'
                        % (p_id, len(pattern), len(p_id2efm_ids[p_id]), ' '.join(sorted(pattern))))
            f.write('\n%d EFMs out of %d were not classified.' % (len(outlier_efm_ids), len(efm_id2efm)))
        else:
            f.write('Failed to find any motifs.\n')

        logging.info('Found %d motifs; %d EFMs out of %d were not classified.'
                     % (len(id2pattern), len(outlier_efm_ids), len(efm_id2efm)))
    logging.info('Saved statistics to %s.' % motif_file)




