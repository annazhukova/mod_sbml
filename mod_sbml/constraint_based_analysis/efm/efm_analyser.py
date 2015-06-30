import logging
import os

import libsbml

from mod_sbml.constraint_based_analysis.efm.acom_classification import acom_classification
from mod_sbml.constraint_based_analysis.efm.control_effective_flux_calculator import classify_efm_by_efficiency
from mod_sbml.constraint_based_analysis.efm.efm_classification import classify_efms
from mod_sbml.constraint_based_analysis.efm.efm_manager import compute_efms
from mod_sbml.constraint_based_analysis.efm.efm_serialization_manager import efm2sbml, serialize_efms_txt, \
    serialize_important_reactions, r_ids2sbml, get_pattern_sorter, serialize_patterns, read_efms, \
    serialize_n_most_effective_efms_txt, pattern2sbml
from mod_sbml.constraint_based_analysis.efm.reaction_classification_by_efm import classify_reactions_by_efm
from mod_sbml.gibbs.reaction_boundary_manager import get_reversible
from mod_sbml.sbml.sbml_manager import reverse_reaction
from mod_sbml.utils.path_manager import create_dirs

__author__ = 'anna'


ZERO_THRESHOLD = 1e-6

def perform_efma(target_r_id, target_r_reversed, sbml, directory, r_id2rev=None, r_ids=None,
                 tree_efm_path="/home/anna/Applications/TreeEFM/tool/TreeEFMseq", max_efm_number=1000, efms=None,
                 threshold=ZERO_THRESHOLD, output_efm_file=None, convert_efms2sbml=False,
                 acom_path='/home/anna/Applications/acom-c/acom-c', min_acom_pattern_len=0, similarity_threshold=0,
                 calculate_patterns=True, min_pattern_len=0, max_pattern_number=10, convert_patterns2sbml=True,
                 calculate_important_reactions=True, imp_rn_threshold=None, rewrite=True,
                 process_sbml=lambda sbml_file, path, header: True, get_file_path=lambda f: f, essential_rn_number=0):
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
                                                      tree_efm_path, r_id2rev, threshold=threshold,
                                                      r_ids=r_ids, rewrite=rewrite)
                if not efms:
                    logging.info('Found no EFMs of interest.')
                    return None, None
            else:
                raise ValueError('No EFMs found :(. Probably, you forgot to specify TreeEFM path?')

        all_id2efm = dict(zip(xrange(1, len(efms) + 1), efms))

        if output_efm_file:
            serialize_efms_txt(all_id2efm, output_efm_file)

    if convert_efms2sbml:
        efm_dir = os.path.join(directory, 'efm_sbml')
        if rewrite or not os.path.exists(efm_dir):
            create_dirs(efm_dir)
            efm2sbml(all_id2efm, directory=efm_dir, sbml=sbml, process_sbml=process_sbml, get_file_path=get_file_path,
                     all_efms_file=output_efm_file, limit=3)

    important_r_ids = None
    if calculate_important_reactions:
        important_r_ids = analyse_important_reactions(all_id2efm, directory, essential_rn_number, get_file_path,
                                                      imp_rn_threshold, model, output_efm_file,
                                                      process_sbml, sbml)

    efficient_efm_ids = analyse_effective_efms(all_id2efm, directory, target_r_id)
    if len(all_id2efm) >= 1000:
        id2efm = {efm_id: efm for (efm_id, efm) in all_id2efm.iteritems() if efm_id in efficient_efm_ids
                  or important_r_ids and not set(efm.to_r_id2coeff().keys()) - important_r_ids}
    else:
        id2efm = all_id2efm

    efm_num = len(id2efm)

    avg_efm_len = sum((len(efm) for efm in id2efm.itervalues())) / efm_num
    min_efm_len = min((len(efm) for efm in id2efm.itervalues()))

    if acom_path:
        analyse_acom(acom_path, avg_efm_len, directory, essential_rn_number, id2efm, min_acom_pattern_len, min_efm_len,
                     r_ids, similarity_threshold)

    if calculate_patterns:
        analyse_patterns(id2efm, avg_efm_len, convert_patterns2sbml, directory, essential_rn_number, get_file_path,
                         max_pattern_number, min_efm_len, min_pattern_len, process_sbml, sbml)

    return id2efm, important_r_ids


def analyse_patterns(id2efm, avg_efm_len, convert_patterns2sbml, directory, essential_rn_number, get_file_path,
                     max_pattern_number, min_efm_len, min_pattern_len, process_sbml, sbml):
    if not min_pattern_len:
        min_pattern_len = min(avg_efm_len / 4, min_efm_len)
        if essential_rn_number:
            min_pattern_len = max(essential_rn_number + 3, min_pattern_len)
    pattern_dir = os.path.join(directory, 'patterns')
    create_dirs(pattern_dir)
    p_id2efm_ids, id2pattern = classify_efms(id2efm, min_pattern_len=min_pattern_len,
                                             max_pattern_num=max_pattern_number)
    sorter = get_pattern_sorter(id2pattern, p_id2efm_ids)
    patterns_txt = os.path.join(pattern_dir, 'patterns.txt')
    serialize_patterns(p_id2efm_ids, id2pattern, patterns_txt, min_pattern_len, sorter=sorter)
    if convert_patterns2sbml:
        pattern_sbml_dir = os.path.join(pattern_dir, 'sbml')
        create_dirs(pattern_sbml_dir)
        pattern2sbml(id2pattern, p_id2efm_ids, pattern_sbml_dir, sbml=sbml, process_sbml=process_sbml,
                     limit=3, get_file_path=get_file_path, all_patterns_file=patterns_txt, total_efm_num=len(id2efm))


def analyse_acom(acom_path, avg_efm_len, directory, essential_rn_number, id2efm, min_acom_pattern_len, min_efm_len,
                 r_ids, similarity_threshold):
    if not min_acom_pattern_len:
        min_acom_pattern_len = min(avg_efm_len / 8, min_efm_len, 12)
        if essential_rn_number:
            min_acom_pattern_len = max(essential_rn_number + 2, min_acom_pattern_len)
    if not similarity_threshold:
        similarity_threshold = int(min_acom_pattern_len * 1.8)
    acom_dir = os.path.join(directory, 'acom')
    create_dirs(acom_dir)
    acom_classification(id2efm, r_ids, directory=acom_dir, acom_path=acom_path,
                        similarity_threshold=similarity_threshold, min_pattern_length=min_acom_pattern_len)


def analyse_effective_efms(all_id2efm, directory, target_r_id):
    id2efficiency = classify_efm_by_efficiency(all_id2efm, target_r_id)
    n = len(id2efficiency) / 3
    efficient_efm_ids = set(sorted(id2efficiency.iterkeys(), key=lambda efm_id: -abs(id2efficiency[efm_id]))[0:n])
    serialize_n_most_effective_efms_txt(all_id2efm, id2efficiency, n=n,
                                        path=os.path.join(directory, 'efficient_efms.txt'))
    return efficient_efm_ids


def analyse_important_reactions(all_id2efm, directory, essential_rn_number, get_file_path, imp_rn_threshold,
                                model, output_efm_file, process_sbml, sbml):
    if imp_rn_threshold is None:
        imp_rn_threshold = len(all_id2efm) / 3
        if essential_rn_number:
            at_least_num = min(int(essential_rn_number * 1.5), len(all_id2efm))
            if at_least_num == len(all_id2efm):
                at_least_num = min(max(len(all_id2efm) / 2, essential_rn_number + 2), len(all_id2efm) - 1)
            imp_rn_threshold = max(at_least_num, imp_rn_threshold)
    rn_dir = os.path.join(directory, 'important')
    create_dirs(rn_dir)
    r_id2efm_ids = classify_reactions_by_efm(all_id2efm)
    if imp_rn_threshold:
        r_id2efm_ids = {r_id: efm_ids for (r_id, efm_ids) in r_id2efm_ids.iteritems()
                        if len(efm_ids) > imp_rn_threshold}
    imp_rn_txt = os.path.join(rn_dir, 'r_list.txt')
    important_r_ids = \
        serialize_important_reactions(r_id2efm_ids, model, imp_rn_txt, imp_rn_threshold)

    def r_updater(r):
        rev_r_id = '-%s' % r.id
        rev = (rev_r_id in r_id2efm_ids and
               (r.id not in r_id2efm_ids or len(r_id2efm_ids[r.id]) < len(r_id2efm_ids[rev_r_id])))
        if rev:
            reverse_reaction(r)
            r.setName(rev_r_id)
        else:
            r.setName(r.id)

    imp_rn_sbml = os.path.join(rn_dir, 'Model_important.xml')
    r_ids2sbml(important_r_ids, sbml, imp_rn_sbml, suffix='important', r_updater=r_updater)
    efm_file_string = '<p>File describing all <b>%d</b> found EFMs is <a href="%s" target="_blank">here</a>.</p>' \
                      % (len(all_id2efm), get_file_path(output_efm_file)) if output_efm_file else ''
    description = \
        '''<h3>Elementary Flux Mode Analysis (EFMA): important reactions</h3>
           <p><b>%d</b> reactions participating in at least <b>%d</b> (out of <b>%d</b>) EFMs are visualised below
           (<a href="%s" target="_blank">detailed description</a>, <a href="%s" download>SBML submodel</a>).</p>
           %s
        ''' % (len(important_r_ids), imp_rn_threshold, len(all_id2efm),
                   get_file_path(imp_rn_txt), get_file_path(imp_rn_sbml), efm_file_string)
    process_sbml(imp_rn_sbml, 'efma_important_reactions', description)
    return important_r_ids
