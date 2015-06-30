import logging
import os
import shutil

from cobra.io.sbml import create_cobra_model_from_sbml_file
import libsbml
from mod_sbml.gibbs.reaction_boundary_manager import get_bounds, set_bounds
from mod_sbml.sbml.sbml_manager import get_products, get_reactants, reverse_reaction
from sbml_vis.graph.graph_properties import ALL_COMPARTMENTS
from mimoza.mimoza_path import JS_SCRIPTS, CSS_SCRIPTS, MIMOZA_FAVICON
from sbml_vis.html.html_generator import create_multi_html
from mimoza_pipeline import process_sbml, get_lib

from mod_sbml.constraint_based_analysis.cobra_constraint_based_analysis.fba_analyser import analyse_by_fba
from mod_sbml.constraint_based_analysis.cobra_constraint_based_analysis.fva_analyser import analyse_by_fva
from mod_sbml.constraint_based_analysis.efm.efm_analyser import perform_efma
from mod_sbml.utils.path_manager import create_dirs
from serialization_manager import get_sbml_r_formula, serialize_model_info
from mod_sbml.sbml.ubiquitous_manager import get_cofactor_m_ids

ZERO_THRESHOLD = 1e-6

__author__ = 'anna'


def constraint_exchange_reactions(model, allowed_exchange_r_id2rev, cofactors=None, min_flux=0.01):
    if not cofactors:
        cofactors = get_cofactor_m_ids(model)

    for r in model.getListOfReactions():
        if r.id in allowed_exchange_r_id2rev:
            rev = allowed_exchange_r_id2rev[r.id]
            if rev:
                # reverse the reaction and set positive bounds,
                # as if both bounds are negative, the glp solver produces an error
                l_b, u_b = get_bounds(r)
                reverse_reaction(r)
                set_bounds(r, -u_b, -l_b)
                allowed_exchange_r_id2rev[r.id] = not rev
            r_l, r_u = get_bounds(r)
            set_bounds(r, max(min_flux, r_l), max(min_flux, r_u))
            continue
        rs, ps = set(get_reactants(r)), set(get_products(r))

        boundary_s_ids = {s_id for s_id in rs if model.getSpecies(s_id).getBoundaryCondition()}
        if boundary_s_ids or not rs:
            r_l, r_u = get_bounds(r)
            # if it's not only cofactors, constrain it
            set_bounds(r, min(r_l, 0), 0 if (ps - cofactors) else r_u)
            continue
        boundary_s_ids = {s_id for s_id in ps if model.getSpecies(s_id).getBoundaryCondition()}
        if boundary_s_ids or not ps:
            r_l, r_u = get_bounds(r)
            # if it's not only cofactors, constrain it
            set_bounds(r, 0 if (rs - cofactors) else r_l, max(r_u, 0))

def analyse_model(sbml, out_r_id, out_rev, res_dir, in_r_id2rev=None, threshold=ZERO_THRESHOLD,
                  do_fva=True, do_fba=True, do_efm=True,
                  efms=None, save_efm_sbml=True, max_efm_number=1000, min_pattern_len=0,
                  max_pattern_num=None, save_pattern_sbml=True, acom_path='/home/anna/Applications/acom-c/acom-c',
                  min_acom_pattern_len=None, imp_rn_threshold=None, rewrite=True, visualise_sbml=True, ub_ch_ids=None,
                  generalize=False, title='', r_ids_of_interest=None, cofactors=None):

    logging.info("Preparing directories...")
    # create directories to store results
    create_dirs(res_dir, False)

    if in_r_id2rev:
        logging.info("Constraining input reactions...")
        doc = libsbml.SBMLReader().readSBML(sbml)
        model = doc.getModel()
        constraint_exchange_reactions(model, allowed_exchange_r_id2rev=in_r_id2rev, cofactors=cofactors)
        sbml = os.path.join(res_dir, '%s_constrained.xml' % os.path.splitext(os.path.basename(sbml))[0])
        libsbml.SBMLWriter().writeSBMLToFile(doc, sbml)

    # copy our model in the result directory
    if os.path.normpath(res_dir) != os.path.normpath(os.path.dirname(sbml)):
        shutil.copy(sbml, res_dir)
        sbml = os.path.join(res_dir, os.path.basename(sbml))

    logging.info("Serializing model info...")
    info_xlsx = os.path.join(res_dir, 'model_info.xlsx')
    serialize_model_info(sbml, info_xlsx)
    in_sbml = sbml
    vis_dir = os.path.join(res_dir, 'visualisation')
    create_dirs(vis_dir, rewrite)
    doc = libsbml.SBMLReader().readSBML(in_sbml)
    model = doc.getModel()
    r_string = lambda r, rev: '<b>%s</b>%s: %s' % (r.getId(), ' (reversed)' if rev else '',
                                                   get_sbml_r_formula(model, r, id=False))
    r = model.getReaction(out_r_id)
    description = '''
        <p>The analysed model is %s (<a href="%s" download>SBML</a>, <a href="%s" download>XLSX</a>).</p>
        <p>Objective reaction is %s.</p>
        %s
    ''' % (model.getName() if model.getName() else model.getId(),
           os.path.join('..', os.path.relpath(in_sbml, res_dir)), os.path.join('..', os.path.relpath(info_xlsx, res_dir)),
           r_string(r, out_rev),
           ('<p>Input reaction%s %s.</p>'
            % (' is' if len(in_r_id2rev) == 1 else 's are',
               '; '.join(r_string(model.getReaction(r_id), rev)
                         for (r_id, rev) in in_r_id2rev.iteritems()))) if in_r_id2rev else '')

    model_data = []

    def vis_sbml(sbml_file, path, description, add=True):
        if visualise_sbml:
            c_id2json_files, c_id2json_vars, c_id2out_c_id, m_id, vis_path, groups_sbml = \
                process_sbml(sbml_file, True, ub_ch_ids=ub_ch_ids, path=path, generalize=generalize)
            if ALL_COMPARTMENTS in c_id2json_vars:
                c_id2json_vars = {ALL_COMPARTMENTS: c_id2json_vars[ALL_COMPARTMENTS]}
                c_id2json_files = {ALL_COMPARTMENTS: c_id2json_files[ALL_COMPARTMENTS]}
            json_files = []
            for f in reduce(lambda l1, l2: l1 + l2, c_id2json_files.itervalues(), []):
                abs_path = os.path.join(vis_path, f)
                common_prefix = os.path.commonprefix([abs_path, res_dir])
                relative_path = os.path.join('..', os.path.relpath(abs_path, common_prefix))
                json_files.append(relative_path)
            if add:
                model_data.append((path, sbml_file, json_files, c_id2json_vars, c_id2out_c_id, m_id, description))
            return groups_sbml
        else:
            return sbml_file

    ess_rn_num = 0
    if do_fva:
        logging.info("Performing FVA...")
        fva_dir = os.path.join(res_dir, 'fva')
        create_dirs(fva_dir, rewrite)

        cobra_model = create_cobra_model_from_sbml_file(sbml)
        r_id2bounds, sbml, opt_val = analyse_by_fva(cobra_model, out_r_id, fva_dir,
                                                    'minimize' if out_rev else 'maximize',
                                                    threshold=threshold, sbml=sbml, rewrite=rewrite,
                                                    process_sbml=vis_sbml, r_ids=r_ids_of_interest,
                                                    get_file_path=\
                                                        lambda f: os.path.join('..', os.path.relpath(f, res_dir)))
        ess_rn_num = len([1 for (r_id, (l, u)) in r_id2bounds.iteritems() if l * u > 0])

    if do_fba:
        logging.info("Performing FBA...")
        fba_dir = os.path.join(res_dir, 'fba')
        create_dirs(fba_dir, rewrite)

        cobra_model = create_cobra_model_from_sbml_file(sbml)
        analyse_by_fba(cobra_model, directory=fba_dir, bm_r_id=out_r_id,
                       objective_sense='minimize' if out_rev else 'maximize', threshold=threshold, sbml=sbml,
                       rewrite=rewrite, process_sbml=vis_sbml, r_ids=r_ids_of_interest,
                       get_file_path=lambda f: os.path.join('..', os.path.relpath(f, res_dir)))

    if do_efm:
        logging.info("Performing EFMA...")
        efm_dir = os.path.join(res_dir, 'efms/')
        create_dirs(efm_dir, rewrite)

        perform_efma(out_r_id, out_rev, sbml=sbml, directory=efm_dir, r_id2rev=in_r_id2rev,
                     max_efm_number=max_efm_number, efms=efms, threshold=threshold,
                     output_efm_file=os.path.join(efm_dir, 'efms.txt'), convert_efms2sbml=save_efm_sbml,
                     acom_path=acom_path, min_acom_pattern_len=min_acom_pattern_len, min_pattern_len=min_pattern_len,
                     max_pattern_number=max_pattern_num, convert_patterns2sbml=save_pattern_sbml,
                     imp_rn_threshold=imp_rn_threshold, rewrite=rewrite, process_sbml=vis_sbml,
                     get_file_path=lambda f: os.path.join('..', os.path.relpath(f, res_dir)),
                     essential_rn_number=ess_rn_num)

    lib_path = os.path.join(vis_dir, 'lib')
    if not os.path.exists(lib_path):
        shutil.copytree(get_lib(), lib_path)

    logging.info('Putting everything together...')

    if not model_data:
        description += '<p>No pathways found.</p>'
    create_multi_html(model_data, title=title, description=description, directory=vis_dir, scripts=JS_SCRIPTS,
                      css=CSS_SCRIPTS, fav=MIMOZA_FAVICON)
