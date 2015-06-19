import logging
import os

import libsbml

from constraint_based_analysis.efm.acom_classification import acom_classification
from constraint_based_analysis.efm.control_effective_flux_calculator import classify_efm_by_efficiency
from constraint_based_analysis.efm.efm_classification import classify_efms
from constraint_based_analysis.efm.efm_manager import compute_efms
from constraint_based_analysis.efm.efm_serialization_manager import efm2sbml, serialize_efms_txt, serialize_efms_xslx, \
    serialize_important_reactions, r_ids2sbml, get_pattern_sorter, serialize_patterns, read_efms, \
    serialize_n_most_effective_efms_txt
from constraint_based_analysis.efm.reaction_classification_by_efm import classify_reactions_by_efm
from gibbs.reaction_boundary_manager import get_bounds, get_reversible
from sbml.sbml_manager import reverse_reaction
from utils.path_manager import create_dirs

__author__ = 'anna'


ZERO_THRESHOLD = 1e-6

def perform_efma(target_r_id, target_r_reversed, r_id2rev_2threshold, sbml, directory, r_ids=None,
                 tree_efm_path="/home/anna/Applications/TreeEFM/tool/TreeEFMseq", max_efm_number=1000, efms=None,
                 threshold=ZERO_THRESHOLD, output_efm_file=None, convert_efms2sbml=False,
                 acom_path='/home/anna/Applications/acom-c/acom-c', min_acom_pattern_len=0, similarity_threshold=0,
                 calculate_patterns=True, min_pattern_len=0, min_efm_num_per_pattern=0, convert_patterns2sbml=True,
                 calculate_important_reactions=True, imp_rn_threshold=None, rewrite=True):
    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()
    if not r_ids:
        r_ids = sorted(r.id for r in model.getListOfReactions())
    rev_r_ids = {r.id for r in (model.getReaction(r_id) for r_id in r_ids) if get_reversible(r)}

    if not rewrite and os.path.exists(output_efm_file) and not efms:
        all_id2efm = read_efms(output_efm_file, r_ids, rev_r_ids)
    else:
        if not efms:
            if tree_efm_path:
                efms, r_ids, rev_r_ids = compute_efms(sbml, directory, max_efm_number, target_r_id, target_r_reversed,
                                                      tree_efm_path, r_id2rev_2threshold, threshold=threshold,
                                                      r_ids=r_ids, rewrite=rewrite)
                if not efms:
                    logging.info('Found no EFMs of interest.')
                    return
            else:
                raise ValueError('No EFMs found :(. Probably, you forgot to specify TreeEFM path?')

        all_id2efm = dict(zip(xrange(1, len(efms) + 1), efms))

        if output_efm_file:
            serialize_efms_txt(all_id2efm, output_efm_file)

    if convert_efms2sbml:
        efm_dir = os.path.join(directory, 'efm_sbml/')
        if rewrite or not os.path.exists(efm_dir):
            create_dirs(efm_dir)
            efm2sbml(all_id2efm, directory=efm_dir, get_name_suffix=lambda efm_id: 'of_len_%d' % len(id2efm[efm_id]),
                     name_prefix='EFM', sbml=sbml)

    important_r_ids = None
    if calculate_important_reactions:
        if imp_rn_threshold is None:
            imp_rn_threshold = len(all_id2efm) / 3

        rn_dir = os.path.join(directory, 'important_reactions/')
        create_dirs(rn_dir)

        r_id2efm_ids = classify_reactions_by_efm(all_id2efm)
        if imp_rn_threshold:
            r_id2efm_ids = {r_id: efm_ids for (r_id, efm_ids) in r_id2efm_ids.iteritems()
                            if len(efm_ids) > imp_rn_threshold}
        important_r_ids = \
            serialize_important_reactions(r_id2efm_ids, model, os.path.join(rn_dir, 'r_list.txt'), imp_rn_threshold)

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

    id2efficiency = classify_efm_by_efficiency(all_id2efm, target_r_id)
    n = min(1000, len(id2efficiency))
    efficient_efm_ids = set(sorted(id2efficiency.iterkeys(), key=lambda efm_id: -abs(id2efficiency[efm_id]))[0:n])
    serialize_n_most_effective_efms_txt(all_id2efm, id2efficiency, n=n,
                                        path=os.path.join(directory, 'efficient_efms.txt'))
    id2efm = {efm_id: efm for (efm_id, efm) in all_id2efm.iteritems() if efm_id in efficient_efm_ids
              or important_r_ids and not set(efm.to_r_id2coeff().keys()) - important_r_ids}

    efm_num = len(id2efm)

    avg_efm_len = sum((len(efm) for efm in id2efm.itervalues())) / efm_num
    min_efm_len = min((len(efm) for efm in id2efm.itervalues()))

    if acom_path:
        if not min_acom_pattern_len:
            min_acom_pattern_len = 2 #min(avg_efm_len / 8, min_efm_len, 12)
        if not similarity_threshold:
            similarity_threshold = int(min_acom_pattern_len * 1.8)

        acom_dir = os.path.join(directory, 'acom/')
        create_dirs(acom_dir)

        acom_classification(id2efm, r_ids, directory=acom_dir, acom_path=acom_path,
                            similarity_threshold=similarity_threshold, min_pattern_length=min_acom_pattern_len)

    if not min_pattern_len:
        min_pattern_len = min(avg_efm_len / 4, min_efm_len)
    if not min_efm_num_per_pattern:
        min_efm_num_per_pattern = efm_num / 6

    if calculate_patterns:
        pattern_dir = os.path.join(directory, 'patterns/')
        create_dirs(pattern_dir)

        p_id2efm_ids, id2pattern = classify_efms(id2efm, min_pattern_len=min_pattern_len,
                                                 min_efm_num=min_efm_num_per_pattern)
        sorter = get_pattern_sorter(id2pattern, p_id2efm_ids, min_pattern_len, min_efm_num_per_pattern)
        serialize_patterns(p_id2efm_ids, id2pattern, os.path.join(pattern_dir, 'patterns.txt'), min_pattern_len,
                           min_efm_num_per_pattern, sorter=sorter)

        if convert_patterns2sbml:
            pattern_sbml_dir = os.path.join(pattern_dir, 'sbml/')
            create_dirs(pattern_sbml_dir)
            efm2sbml(id2pattern, directory=pattern_sbml_dir,
                     get_name_suffix=\
                         lambda p_id: 'of_%d_EFMs_of_len_%d' % (len(p_id2efm_ids[p_id]), len(id2pattern[p_id])),
                     name_prefix='Pattern', sbml=sbml)

    return all_id2efm, important_r_ids
