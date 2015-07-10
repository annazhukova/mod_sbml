from collections import defaultdict
import logging
import os
import shutil

from cobra.io.sbml import create_cobra_model_from_sbml_file
import libsbml

from cobra_constraint_based_analysis.fba_manager import serialize_fva, serialize_fluxes
from cobra_constraint_based_analysis.model_manager import format_r_id
from efm.clique_detection import detect_cliques
from efm.efm_classification import classify_efms
from efm.efm_serialization_manager import r_ids2sbml, serialize_efms_txt, serialize_efm, \
    serialize_important_reactions, get_pattern_sorter, serialize_patterns, serialize_pattern, \
    serialize_cliques, serialize_clique
from efm.reaction_classification_by_efm import get_important_reactions
from mod_sbml.gibbs.reaction_boundary_manager import get_bounds, set_bounds
from mod_sbml.sbml.sbml_manager import get_products, get_reactants, reverse_reaction, get_r_comps
from mimoza_pipeline import process_sbml
from mod_sbml.constraint_based_analysis.cobra_constraint_based_analysis.fba_analyser import analyse_by_fba, \
    create_fba_model
from mod_sbml.constraint_based_analysis.cobra_constraint_based_analysis.fva_analyser import analyse_by_fva, \
    create_fva_model, create_essential_r_ids_model
from mod_sbml.constraint_based_analysis.efm.efm_analyser import calculate_imp_rn_threshold, get_efms, \
    calculate_min_pattern_len, calculate_min_clique_len
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

def analyse_model(sbml, out_r_id, out_rev, res_dir, in_r_id2rev=None, threshold=ZERO_THRESHOLD, do_fva=True,
                  do_fba=True, do_efm=True, efms=None, max_efm_number=1000, min_pattern_len=0,
                  imp_rn_threshold=None, ub_ch_ids=None,
                  generalize=False, title='', r_ids_of_interest=None, cofactors=None):
    logging.info("Preparing directories...")
    # create directories to store results
    create_dirs(res_dir, False)
    get_f_path = lambda f: os.path.join('..', os.path.relpath(f, res_dir))

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
    doc = libsbml.SBMLReader().readSBML(in_sbml)
    model = doc.getModel()
    r_string = lambda r, rev: '<b>%s</b>%s: %s' % (r.getId(), ' (reversed)' if rev else '',
                                                   get_sbml_r_formula(model, r, id=False))
    r = model.getReaction(out_r_id)
    description = '''
        <div class="mediabox">
            <h3>Input data</h3>
            <p>The analysed model is %s (<a href="%s" download>SBML</a>, <a href="%s" download>XLSX</a>).</p>
            <p>Objective reaction is %s.</p>
            %s
        </div>
    ''' % (model.getName() if model.getName() else model.getId(),
           os.path.join('..', os.path.relpath(in_sbml, res_dir)), os.path.join('..', os.path.relpath(info_xlsx, res_dir)),
           r_string(r, out_rev),
           ('<p>Input reaction%s %s.</p>'
            % (' is' if len(in_r_id2rev) == 1 else 's are',
               '; '.join(r_string(model.getReaction(r_id), rev)
                         for (r_id, rev) in in_r_id2rev.iteritems()))) if in_r_id2rev else '')

    id2mask = defaultdict(lambda: 0)
    layer2mask = {}
    mask_shift = 4

    ess_rn_num = 0
    objective_sense = 'minimize' if out_rev else 'maximize'
    vis_r_ids = set()
    if do_fva:
        logging.info("Performing FVA...")
        fva_dir = os.path.join(res_dir, 'fva')
        create_dirs(fva_dir)

        cobra_model = create_cobra_model_from_sbml_file(sbml)
        r_id2bounds, opt_val = analyse_by_fva(cobra_model=cobra_model, bm_r_id=out_r_id,
                                              objective_sense=objective_sense, threshold=threshold)
        if opt_val:
            fva_file = os.path.join(fva_dir, 'fva.txt')
            caption = '%s reaction %s (%s): %.4g\n==================================\n'\
                      % (objective_sense, out_r_id,
                         get_sbml_r_formula(model, model.getReaction(out_r_id), id=True, comp=False), opt_val)
            ess_rn_num, var_rn_num = serialize_fva(cobra_model, r_id2bounds, fva_file, r_ids=r_ids_of_interest,
                                                   title=caption)
            essential_r_id2rev = {r_id: u < 0 for (r_id, (l, u)) in r_id2bounds.iteritems() if l * u > 0}
            ess_r_ids = set(essential_r_id2rev.keys())
            update_vis_layers(ess_r_ids, 'FVA: essential', id2mask, layer2mask, mask_shift, model, vis_r_ids)
            mask_shift += 1

            if opt_val and sbml:
                efm_sbml = os.path.join(fva_dir, 'Model_FVA.xml')
                sbml = create_fva_model(sbml, r_id2bounds, efm_sbml)
                ess_sbml = create_essential_r_ids_model(sbml, essential_r_id2rev, fva_dir)
                description += \
                    '''<div class="mediabox">
                           <h3>Flux Variability Analysis (FVA)</h3>
                           <p>Optimal flux through the objective reaction is <b>%g</b>.</p>
                           <p>It requires <b>%d</b> essential and <b>%d</b> variable reactions
                           (<a href="%s" target="_blank">detailed description</a>, <a href="%s" download>SBML submodel of essential reactions</a>,
                           <a href="%s" download>SBML submodel of essential and variable reactions</a>).</p>
                       </div>
                    ''' % (opt_val, ess_rn_num, var_rn_num, get_f_path(fva_file), get_f_path(ess_sbml), get_f_path(efm_sbml))

    if do_fba:
        logging.info("Performing FBA...")
        fba_dir = os.path.join(res_dir, 'fba')
        create_dirs(fba_dir)

        cobra_model = create_cobra_model_from_sbml_file(sbml)
        r_id2val, opt_val = analyse_by_fba(cobra_model, bm_r_id=out_r_id, objective_sense=objective_sense,
                                           threshold=threshold)
        if opt_val:
            fba_r_ids = set(r_id2val.keys())
            update_vis_layers(fba_r_ids, 'FBA', id2mask, layer2mask, mask_shift, model, vis_r_ids)
            mask_shift += 1

            fba_file = os.path.join(fba_dir, 'fba.txt')
            caption = '%s reaction %s (%s): %.4g\n==================================\n'\
                      % (objective_sense, out_r_id,
                         get_sbml_r_formula(model, model.getReaction(out_r_id), id=True, comp=False), opt_val)
            serialize_fluxes(cobra_model, r_id2val, path=fba_file, r_ids=r_ids_of_interest, title=caption)

            if opt_val and sbml:
                efm_sbml = os.path.join(fba_dir, 'Model_FBA.xml')
                create_fba_model(sbml, r_id2val, efm_sbml)
                description += \
                    '''<div class="mediabox">
                           <h3>Flux Balance Analysis (FBA)</h3>
                           <p>Optimal flux through the objective reaction is <b>%g</b>.</p>
                           <p>It involves <b>%d</b> reactions
                           (<a href="%s" target="_blank">detailed description</a>, <a href="%s" download>SBML submodel</a>).</p>
                        </div>
                    ''' % (opt_val, len(r_id2val), get_f_path(fba_file), get_f_path(efm_sbml))

    if do_efm:
        logging.info("Performing EFMA...")
        efm_dir = os.path.join(res_dir, 'efms/')
        create_dirs(efm_dir)

        id2efm = get_efms(target_r_id=out_r_id, target_r_reversed=out_rev, r_id2rev=in_r_id2rev, sbml=sbml,
                          directory=efm_dir, max_efm_number=max_efm_number, threshold=threshold, efms=efms,
                          rewrite=out_r_id, r_ids=r_ids_of_interest)
        if id2efm:
            all_efm_intersection = reduce(lambda p1, p2: p1.intersection(p2), id2efm.itervalues(),
                                          next(id2efm.itervalues()))
            all_efm_file = os.path.join(efm_dir, 'efms.txt')
            serialize_efms_txt(id2efm, all_efm_file, out_r_id, all_efm_intersection)
            description += \
                '''<div class="mediabox">
                       <h3>Elementary Flux Mode Analysis (EFMA)</h3>
                       <p>Calculated <b>%d</b> EFMs (<a href="%s" target="_blank">detailed description</a>).</p>
                ''' % (len(id2efm), get_f_path(all_efm_file))

            # 3 shortest EFMs
            mask_shift, description\
                = process_n_fms(sbml, model, description, directory=efm_dir, id2fm=id2efm, id2mask=id2mask,
                                       layer2mask=layer2mask, mask_shift=mask_shift, vis_r_ids=vis_r_ids,
                                       sort_criterion='shortest', name='EFM', sorter=lambda (_, efm): len(efm),
                                       serializer=lambda efm_id, efm_txt:
                                       serialize_efm(efm_id, id2efm[efm_id], model, efm_txt, out_r_id),
                                       suffix=lambda _: '', get_f_path=get_f_path)

            # Important reactions
            imp_rn_threshold = calculate_imp_rn_threshold(id2efm, ess_rn_num, imp_rn_threshold)
            rn_dir = os.path.join(efm_dir, 'important')
            create_dirs(rn_dir)
            r_id2efm_ids, important_r_ids = get_important_reactions(id2efm, imp_rn_threshold)
            imp_rn_txt = os.path.join(rn_dir, 'r_list.txt')
            serialize_important_reactions(r_id2efm_ids, model, imp_rn_txt, imp_rn_threshold)
            imp_rn_sbml = os.path.join(rn_dir, 'Model_important.xml')
            r_ids2sbml(important_r_ids, sbml, imp_rn_sbml, suffix='important')

            description += '''<p>Found <b>%d</b> reactions participating in at least <b>%d</b> EFMs
                   (<a href="%s" target="_blank">detailed description</a>, <a href="%s" download>SBML submodel</a>).</p>
                   ''' % (len(important_r_ids), imp_rn_threshold, get_f_path(imp_rn_txt), get_f_path(imp_rn_sbml))

            avg_efm_len = sum((len(efm) for efm in id2efm.itervalues())) / len(id2efm)
            min_efm_len = min((len(efm) for efm in id2efm.itervalues()))

            description += '</div>'

            description += '''<div class="mediabox">
                       <h3>Patterns</h3>'''

            # Patterns
            min_pattern_len = calculate_min_pattern_len(avg_efm_len, ess_rn_num, min_efm_len, min_pattern_len)
            pattern_dir = os.path.join(efm_dir, 'patterns')
            create_dirs(pattern_dir)
            p_id2efm_ids, id2pattern = classify_efms(id2efm, min_pattern_len=min_pattern_len,
                                                                           min_efm_num=imp_rn_threshold)
            sorter = get_pattern_sorter(id2pattern, p_id2efm_ids)
            patterns_txt = os.path.join(pattern_dir, 'patterns.txt')
            serialize_patterns(p_id2efm_ids, id2pattern, patterns_txt, min_pattern_len, all_efm_intersection,
                               min_efm_num=imp_rn_threshold, sorter=sorter)

            description += '''
                       <p>Detected <b>%d</b> patterns of length at least <b>%d</b>,
            found in at least <b>%d</b> EFMs
            (<a href="%s" target="_blank">detailed description</a>).</p>
            ''' % (len(id2pattern), min_pattern_len, imp_rn_threshold, get_f_path(patterns_txt))

            # 3 most common patterns
            mask_shift, description\
                = process_n_fms(sbml, model, description, directory=pattern_dir, id2fm=id2pattern, id2mask=id2mask,
                                       layer2mask=layer2mask, mask_shift=mask_shift, vis_r_ids=vis_r_ids,
                                       sort_criterion='most common', name='pattern',
                                       sorter=lambda (p_id, _): -len(p_id2efm_ids[p_id]),
                                       serializer=lambda p_id, pattern_txt:
                                       serialize_pattern(p_id, id2pattern[p_id], p_id2efm_ids[p_id], model, pattern_txt),
                                       suffix=lambda p_id:
                                       'found in <b>%d</b> EFMs ' % len(p_id2efm_ids[p_id]),
                                       get_f_path=get_f_path)
            description += '</div>'

            # Cliques
            description += '''<div class="mediabox">
                       <h3>Cliques</h3>'''
            min_clique_len = calculate_min_clique_len(avg_efm_len, ess_rn_num, min_clique_len=None,
                                                      min_efm_len=imp_rn_threshold)
            clique_dir = os.path.join(efm_dir, 'cliques')
            create_dirs(clique_dir)
            id2clique = detect_cliques(id2efm, min_clique_size=min_clique_len, efm_num=imp_rn_threshold)
            cliques_txt = os.path.join(clique_dir, 'cliques.txt')
            serialize_cliques(id2clique, cliques_txt, min_clique_len, all_efm_intersection, min_efm_num=imp_rn_threshold)
            description += '''<p>Detected <b>%d</b> cliques of length at least <b>%d</b>,
            with <b>%d</b> as min number of common EFMs for related reactions
            (<a href="%s" target="_blank">detailed description</a>).</p>''' \
                           % (len(id2clique), min_clique_len, imp_rn_threshold, get_f_path(cliques_txt))

            # 3 longest cliques
            mask_shift, description \
                = process_n_fms(sbml, model, description, directory=clique_dir, id2fm=id2clique, id2mask=id2mask,
                                       layer2mask=layer2mask, mask_shift=mask_shift, vis_r_ids=vis_r_ids,
                                       sort_criterion='longest', name='clique',
                                       sorter=lambda (_, cl): -len(cl),
                                       serializer=lambda cl_id, clique_txt:
                                       serialize_clique(cl_id, id2clique[cl_id], model, clique_txt),
                                       suffix=lambda _: '', get_f_path=get_f_path)

            description += '</div>'

    logging.info('Putting everything together...')

    if not vis_r_ids:
        description += '''<div class="mediabox">
            <h3>No pathways found :(</div>'''

    combined_sbml = os.path.join(res_dir, 'Combined_model.xml')
    r_ids2sbml(vis_r_ids, in_sbml, combined_sbml, 'combined')
    process_sbml(combined_sbml, True, ub_ch_ids=ub_ch_ids, path='visualization', generalize=generalize,
                 id2mask=id2mask, layer2mask=layer2mask, tab2html={'Analysis': description}, title=title)


def update_vis_layers(r_ids, layer, id2mask, layer2mask, mask_shift, model, vis_r_ids):
    r_ids |= {format_r_id(r_id, False) for r_id in r_ids}
    vis_r_ids |= r_ids
    l_mask = 1 << mask_shift
    layer2mask[layer] = l_mask
    for r_id in r_ids:
        r = model.getReaction(r_id)
        if r:
            id2mask[r_id] |= l_mask
            id2mask.update({c_id: id2mask[c_id] | l_mask for c_id in get_r_comps(r_id, model)})
            id2mask.update({s_id: id2mask[s_id] | l_mask for s_id in get_reactants(r)})
            id2mask.update({s_id: id2mask[s_id] | l_mask for s_id in get_products(r)})


def process_n_fms(sbml, model, description, directory, id2fm, id2mask, layer2mask, mask_shift, vis_r_ids,
                  sort_criterion, name, sorter, serializer, suffix, get_f_path):
    limit = min(3, len(id2fm))
    if limit:
        sbml_dir = os.path.join(directory, 'sbml')
        create_dirs(sbml_dir)
        description += '''<p>%d %s %ss are:</p>
        <ol class="small">''' % (limit, sort_criterion, name)
        for fm_id, fm in sorted(id2fm.iteritems(), key=sorter):
            r_id2coeff = fm.to_r_id2coeff()
            c_name = name[0].upper() + name[1:]
            fm_name = '%s_%d_of_len_%d' % (c_name, fm_id, len(fm))
            fm_sbml = os.path.join(sbml_dir, '%s.xml' % fm_name)
            r_ids2sbml(r_id2coeff.keys(), sbml, fm_sbml, fm_name)
            fm_txt = os.path.join(sbml_dir, '%s.txt' % fm_name)
            serializer(fm_id, fm_txt)

            description += '''<li  class="small">%s %d of length <b>%d</b> %s
                   (<a href="%s" target="_blank">detailed description</a>, <a href="%s" download>SBML submodel</a>);</li>
                ''' % (c_name, fm_id, len(fm), suffix(fm_id), get_f_path(fm_txt), get_f_path(fm_sbml))

            update_vis_layers(set(r_id2coeff), '%s %d' % (c_name, fm_id), id2mask, layer2mask, mask_shift, model, vis_r_ids)
            mask_shift += 1

            limit -= 1
            if limit <= 0:
                break
        description += '''</ol>
        </p>'''
    return mask_shift, description

