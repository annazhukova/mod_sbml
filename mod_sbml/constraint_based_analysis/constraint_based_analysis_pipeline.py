import logging
import os
import shutil

from cobra.io.sbml import create_cobra_model_from_sbml_file
from mimoza.mimoza_path import JS_SCRIPTS, CSS_SCRIPTS, MIMOZA_FAVICON
from sbml_vis.html.html_generator import create_multi_html
from mimoza_pipeline import process_sbml, get_lib

from mod_sbml.constraint_based_analysis.cobra_constraint_based_analysis.fba_analyser import analyse_by_fba
from mod_sbml.constraint_based_analysis.cobra_constraint_based_analysis.fva_analyser import analyse_by_fva
from mod_sbml.constraint_based_analysis.efm.efm_analyser import perform_efma
from mod_sbml.utils.path_manager import create_dirs

ZERO_THRESHOLD = 1e-6

__author__ = 'anna'

def analyse_model(sbml, out_r_id, out_rev, in_r_id2rev_2threshold, res_dir, constraint_cobra_reactions=None,
                  threshold=ZERO_THRESHOLD, do_fva=True, do_fba=True, do_efm=True, efms=None, save_efm_sbml=False,
                  max_efm_number=1000, min_pattern_len=0, max_pattern_num=10, save_pattern_sbml=True,
                  min_acom_pattern_len=None, imp_rn_threshold=None, rewrite=True, visualise_sbml=True,
                  ub_ch_ids=None, generalize=False, title=''):
    model_data = []

    def vis_sbml(sbml_file, path):
        if visualise_sbml:
            json_files, c_id2json_vars, c_id2out_c_id, m_id, vis_path, groups_sbml = \
                process_sbml(sbml_file, True, ub_ch_ids=ub_ch_ids, path=path, generalize=generalize)
            json_files = ['file://%s/%s' % (vis_path, f) for f in json_files]
            model_data.append((path, sbml_file, json_files, c_id2json_vars, c_id2out_c_id, m_id))
            return groups_sbml
        else:
            return sbml_file

    logging.info("Preparing directories...")
    # create directories to store results
    create_dirs(res_dir, False)
    # copy our model there
    if os.path.normpath(res_dir) != os.path.normpath(os.path.dirname(sbml)):
        shutil.copy(sbml, res_dir)
        sbml = os.path.join(res_dir, os.path.basename(sbml))

    if do_fva:
        logging.info("Performing FVA...")
        fva_dir = os.path.join(res_dir, 'fva')
        create_dirs(fva_dir, rewrite)

        cobra_model = create_cobra_model_from_sbml_file(sbml)
        if constraint_cobra_reactions:
            constraint_cobra_reactions(cobra_model)

        _, sbml = analyse_by_fva(cobra_model, out_r_id, fva_dir, 'minimize' if out_rev else 'maximize',
                                 threshold=threshold, sbml=sbml, rewrite=rewrite,
                                 process_sbml=vis_sbml)

    if do_fba:
        logging.info("Performing FBA...")
        fba_dir = os.path.join(res_dir, 'fba')
        create_dirs(fba_dir, rewrite)

        cobra_model = create_cobra_model_from_sbml_file(sbml)
        if constraint_cobra_reactions:
            constraint_cobra_reactions(cobra_model)

        # FBA
        analyse_by_fba(cobra_model, directory=fba_dir,
                       bm_r_id=out_r_id, objective_sense='minimize' if out_rev else 'maximize',
                       threshold=threshold, sbml=sbml, rewrite=rewrite,
                       process_sbml=vis_sbml)

    if do_efm:
        logging.info("Performing EFMA...")
        efm_dir = os.path.join(res_dir, 'efms/')
        create_dirs(efm_dir, rewrite)

        perform_efma(out_r_id, out_rev, r_id2rev_2threshold=in_r_id2rev_2threshold, sbml=sbml, directory=efm_dir,
                     max_efm_number=max_efm_number, threshold=threshold, efms=efms,
                     output_efm_file=os.path.join(efm_dir, 'efms.txt'), min_pattern_len=min_pattern_len,
                     max_pattern_number=max_pattern_num, min_acom_pattern_len=min_acom_pattern_len,
                     convert_efms2sbml=save_efm_sbml,
                     convert_patterns2sbml=save_pattern_sbml, imp_rn_threshold=imp_rn_threshold, rewrite=rewrite,
                     process_sbml=vis_sbml)

    vis_dir = os.path.join(res_dir, 'visualisation/')
    create_dirs(vis_dir, rewrite)

    lib_path = os.path.join(vis_dir, 'lib')
    if not os.path.exists(lib_path):
        shutil.copytree(get_lib(), lib_path)
    create_multi_html(model_data, title=title, statistics='Performed COBRA. <br>See results in the corresponding tabs.',
                      directory=vis_dir, scripts=JS_SCRIPTS,
                      css=CSS_SCRIPTS, fav=MIMOZA_FAVICON)
