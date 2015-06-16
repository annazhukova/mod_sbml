import logging
import os

import libsbml

from constraint_based_analysis.efm.acom_classification import acom_classification
from constraint_based_analysis.efm.efm_classification import classify_efms
from constraint_based_analysis.efm.efm_manager import compute_efms, get_binary_efm_length
from constraint_based_analysis.efm.efm_serialization_manager import efm2sbml, serialize_efms_txt, serialize_efms_xslx, \
    serialize_important_reactions, r_ids2sbml, get_pattern_sorter, serialize_patterns
from constraint_based_analysis.efm.reaction_classification_by_efm import classify_reactions_by_efm
from sbml.sbml_manager import reverse_reaction
from utils.path_manager import create_dirs

__author__ = 'anna'


ZERO_THRESHOLD = 1e-6

def perform_efma(in_r_id, in_r_reversed, out_r_id2rev_2threshold, sbml, directory, r_ids=None,
                 tree_efm_path="/home/anna/Applications/TreeEFM/tool/TreeEFMseq", max_efm_number=1000, efms=None,
                 threshold=ZERO_THRESHOLD, output_efm_file=None, convert_efms2sbml=False,
                 acom_path='/home/anna/Applications/acom-c/acom-c', min_acom_pattern_len=0, similarity_threshold=0,
                 calculate_patterns=True, min_pattern_len=0, min_efm_num_per_pattern=0, convert_patterns2sbml=True,
                 calculate_important_reactions=True, imp_rn_threshold=None):
    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()

    if tree_efm_path:
        efms, r_ids, rev_r_ids = compute_efms(sbml, directory, max_efm_number, in_r_id, in_r_reversed, tree_efm_path,
                                              out_r_id2rev_2threshold, threshold=threshold, r_ids=r_ids)
    else:
        if not efms:
            raise ValueError('No EFMs found :(. Probably, you forgot to specify TreeEFM path?')
        rev_r_ids = {r_id for r_id in r_ids if model.getReaction(r_id).getReversible()}

    if not efms:
        logging.info('Found no EFMs of interest.')
        return

    id2efm = dict(zip(xrange(1, len(efms) + 1), efms))

    if output_efm_file:
        if -1 != output_efm_file.find('.xslx'):
            serialize_efms_xslx(sbml, efms, r_ids, rev_r_ids, output_efm_file)
        else:
            serialize_efms_txt(sbml, id2efm, r_ids, rev_r_ids, output_efm_file)

    if convert_efms2sbml:
        efm_dir = os.path.join(directory, 'efm_sbml/')
        create_dirs(efm_dir)
        efm2sbml(id2efm, directory=efm_dir,
                 get_name_suffix=lambda efm_id: 'of_len_%d' % get_binary_efm_length(id2efm[efm_id][0]),
                 name_prefix='EFM', sbml=sbml, r_ids=r_ids,
                 rev_r_ids=rev_r_ids, binary=False)

    avg_efm_len = sum((len(efm[1]) for efm in efms)) / len(efms)

    if acom_path:
        if not min_acom_pattern_len:
            min_acom_pattern_len = avg_efm_len / 4
        if not similarity_threshold:
            similarity_threshold = 2 * min_acom_pattern_len

        acom_dir = os.path.join(directory, 'acom/')
        create_dirs(acom_dir)

        acom_classification(id2efm, r_ids, rev_r_ids, directory=acom_dir, acom_path=acom_path,
                            similarity_threshold=similarity_threshold, min_pattern_length=min_acom_pattern_len)

    if not min_pattern_len:
        min_pattern_len = avg_efm_len / 4
    if not min_efm_num_per_pattern:
        min_efm_num_per_pattern = len(efms) / 5

    if calculate_important_reactions:
        if imp_rn_threshold is None:
            imp_rn_threshold = min_efm_num_per_pattern

        rn_dir = os.path.join(directory, 'important_reactions/')
        create_dirs(rn_dir)

        r_id2efm_ids = classify_reactions_by_efm(id2efm, r_ids, rev_r_ids)
        r_id2efm_ids = {r_id: efm_ids for (r_id, efm_ids) in r_id2efm_ids.iteritems()
                        if len(efm_ids) > imp_rn_threshold}
        important_r_ids = {r_id[1:] if '-' == r_id[0] else r_id for r_id in r_id2efm_ids.iterkeys()}
        serialize_important_reactions(r_id2efm_ids, model, os.path.join(rn_dir, 'r_list.txt'))

        def r_updater(r):
            rev_r_id = '-%s' % r.id
            rev = (rev_r_id in r_id2efm_ids and
                   (r.id not in r_id2efm_ids or len(r_id2efm_ids[r.id]) < len(r_id2efm_ids[rev_r_id])))
            if rev:
                reverse_reaction(r)
                r.setName(rev_r_id)
            else:
                r.setName(r.id)

        r_ids2sbml(important_r_ids, sbml, '%s/Model_important.xml' % rn_dir,
                   suffix='important', r_updater=r_updater)

    if calculate_patterns:
        pattern_dir = os.path.join(directory, 'patterns/')
        create_dirs(pattern_dir)

        p_id2efm_ids, id2pattern = classify_efms(efms, min_pattern_len=min_pattern_len,
                                                 min_efm_num=min_efm_num_per_pattern)
        sorter = get_pattern_sorter(id2pattern, p_id2efm_ids, min_pattern_len, min_efm_num_per_pattern)
        serialize_patterns(p_id2efm_ids, id2pattern, r_ids, rev_r_ids,
                           os.path.join(pattern_dir, 'patterns.txt'), sorter=sorter)

        if convert_patterns2sbml:
            pattern_sbml_dir = os.path.join(pattern_dir, 'sbml/')
            create_dirs(pattern_sbml_dir)
            efm2sbml(id2pattern, directory=pattern_sbml_dir,
                     get_name_suffix=\
                         lambda p_id: 'of_%d_EFMs_of_len_%d'
                                      % (len(p_id2efm_ids[p_id]), get_binary_efm_length(id2pattern[p_id])),
                     name_prefix='Pattern', sbml=sbml, r_ids=r_ids, rev_r_ids=rev_r_ids, binary=True)
