import libsbml
from em.acom_classification import acom_classification
from em.em_classification_binary import classify_efms
from em.em_manager import compute_efms, get_binary_efm_length
from em.em_serialization_manager import basic_r_style, serialize_efms_xslx, serialize_efms_txt, efm2sbml, \
    get_pattern_sorter, serialize_patterns

__author__ = 'anna'


ZERO_THRESHOLD = 1e-9

def perform_efma(sbml, in_r_id, in_r_reversed, out_r_id2rev_2threshold, em_number,
                 min_pattern_len, model_dir, output_pattern_file=None, output_efm_file=None,
                 pattern_dir=None, efm_dir=None, over_expressed_r_ids=set(),
                 under_expressed_r_ids=set(),
                 threshold=ZERO_THRESHOLD, r_ids=None, r_id2style=basic_r_style,
                 tree_efm_path="/home/anna/Applications/TreeEFM/tool/TreeEFMseq", min_efm_num_per_pattern=2):
    efms, r_ids = compute_efms(sbml, model_dir, em_number, in_r_id, in_r_reversed, tree_efm_path,
                               out_r_id2rev_2threshold,
                               threshold=threshold, r_ids=r_ids)

    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()
    reversible_r_ids = {r_id for r_id in r_ids if model.getReaction(r_id).getReversible()}

    if output_efm_file:
        if -1 != output_efm_file.find('.xsls'):
            serialize_efms_xslx(sbml, efms, r_ids, reversible_r_ids, output_efm_file,
                                over_expressed_r_ids=over_expressed_r_ids,
                                under_expressed_r_ids=under_expressed_r_ids, r_id2style=r_id2style)
        else:
            serialize_efms_txt(sbml, efms, r_ids, reversible_r_ids, output_efm_file)

    if efm_dir is not None:
        id2efm = dict(zip(xrange(0, len(efms)), efms))
        efm2sbml(id2efm, directory=efm_dir,
                 get_name_suffix=lambda efm_id: 'of_len_%d' % get_binary_efm_length(id2efm[efm_id][0]),
                 name_prefix='EFM', sbml=sbml, r_ids=r_ids,
                 rev_r_ids=reversible_r_ids, binary=False)

    if efms:
        if not min_pattern_len:
            min_pattern_len = 2 * (sum((len(efm[1]) for efm in efms)) / len(efms)) / 3
        if not min_efm_num_per_pattern:
            min_efm_num_per_pattern = len(efms) / 5

        acom_classification(efms, r_ids, reversible_r_ids, directory=model_dir,
                            similarity_threshold=min_pattern_len + min_pattern_len / 3,
                            min_pattern_size=min_pattern_len, acom_path='/home/anna/Applications/acom-c/acom-c')

        p_id2efm_ids, id2pattern = classify_efms(efms, min_pattern_len=min_pattern_len, min_efm_num=min_efm_num_per_pattern)

        if output_pattern_file:
            sorter = get_pattern_sorter(id2pattern, p_id2efm_ids, min_pattern_len, min_efm_num_per_pattern)
            serialize_patterns(p_id2efm_ids, id2pattern, r_ids, reversible_r_ids, output_pattern_file, sorter=sorter)

        if pattern_dir is not None:
            efm2sbml(id2pattern, directory=pattern_dir,
                     get_name_suffix=lambda p_id: 'of_%d_EFMs_of_len_%d'
                                                  % (len(p_id2efm_ids[p_id]), get_binary_efm_length(id2pattern[p_id])),
                     name_prefix='Pattern', sbml=sbml, r_ids=r_ids,
                     rev_r_ids=reversible_r_ids, binary=True)

