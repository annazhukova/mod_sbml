import logging
import os
from mod_sbml.constraint_based_analysis.cobra_constraint_based_analysis.fba_manager import optimise_biomass, \
    get_fluxes_larger_than_threshold, serialize_fluxes, show_reaction_values, round_value

from mod_sbml.constraint_based_analysis.efm.efm_serialization_manager import r_ids2sbml

from mod_sbml.constraint_based_analysis.cobra_constraint_based_analysis.model_manager import format_r_id

from mod_sbml.sbml.sbml_manager import reverse_reaction
from mod_sbml.serialization.serialization_manager import get_cobra_r_formula

__author__ = 'anna'


def create_fba_model(sbml, r_id2val, new_sbml):
    def r_updater(r):
        val = r_id2val[format_r_id(r.id)]
        if val < 0:
            reverse_reaction(r)
            r.setName('%g -%s' % (-val, r.id))
        else:
            r.setName('%g %s' % (val, r.id))

    r_ids = set(r_id2val.iterkeys()) | {format_r_id(r_id, False) for r_id in r_id2val.iterkeys()}
    r_ids2sbml(r_ids, sbml, new_sbml, 'FBA', r_updater)
    return new_sbml


def analyse_by_fba(cobra_model, bm_r_id, directory, objective_sense='maximize', threshold=0,
                   fig_path=None, var_r_id=None, r_ids=None, sbml=None, rewrite=True,
                   process_sbml=lambda sbml_file, path: True, get_file_path=lambda f: f):
    fba_file = os.path.join(directory, 'fba.txt')
    if not rewrite and os.path.exists(fba_file) and (not fig_path or os.path.exists(fig_path)):
        if sbml:
            new_sbml = os.path.join(directory, 'Model_FBA.xml')
            if os.path.exists(new_sbml):
                return new_sbml
        else:
            return sbml

    cobra_bm_r_id = format_r_id(bm_r_id)
    opt_value = optimise_biomass(cobra_model, cobra_bm_r_id, objective_sense, level=logging.DEBUG)
    opt_value = round_value(opt_value)
    r_id2val = get_fluxes_larger_than_threshold(cobra_model, threshold=threshold)
    title = '%s reaction %s (%s): %.4g\n==================================\n'\
            % (objective_sense, bm_r_id, get_cobra_r_formula(cobra_model.reactions.get_by_id(cobra_bm_r_id)),
               opt_value)
    serialize_fluxes(cobra_model, r_id2val, path=fba_file, r_ids=r_ids, title=title)

    if fig_path:
        var_r_id = format_r_id(var_r_id)
        var_r = cobra_model.reactions.get_by_id(var_r_id)
        show_reaction_values(cobra_model, var_r_id, {bm_r_id}, bm_r_id,
                             fig_path, "%s flux (micromol/min/gDW) bounds" % var_r.name,
                             "%s production (micromol/min/gDW)" % bm_r_id,
                             "Effect of varying %s on %s" % (var_r.name, bm_r_id),
                             var_r.upper_bound * 2.5, minimized=False, constrained=False)
    if opt_value and sbml:
        new_sbml = os.path.join(directory, 'Model_FBA.xml')
        sbml = create_fba_model(sbml, r_id2val, new_sbml)
        description = \
            '''<h3>Flux Balance Analysis (FBA)</h3>
               <p>Optimal flux through the objective reaction is <b>%g</b>.</p>
               <p>It involves <b>%d</b> reactions
               (<a href="%s" target="_blank">detailed description</a>, <a href="%s" download>SBML submodel</a>)
               that are visualised below.</p>
            ''' % (opt_value, len(r_id2val), get_file_path(fba_file), get_file_path(new_sbml))
        process_sbml(sbml, 'fba', description)
    return sbml, opt_value