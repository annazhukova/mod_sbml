import logging
import os
import libsbml
from constraint_based_analysis.cobra_constraint_based_analysis.fba_manager import get_fluxes_larger_than_threshold, \
    optimise_biomass, show_reaction_values, serialize_fluxes
from constraint_based_analysis.cobra_constraint_based_analysis.model_manager import format_r_id
from sbml.sbml_manager import submodel
from serialization.serialization_manager import get_cobra_r_formula

__author__ = 'anna'


def create_fba_model(sbml, r_id2val, res_dir):
    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()
    submodel(set(r_id2val.iterkeys()) | {format_r_id(r_id, False) for r_id in r_id2val.iterkeys()} |
             {'R_' + format_r_id(r_id, False) for r_id in r_id2val.iterkeys()}, model)
    model.setId('%s_FBA' % model.getId())
    model.setName('%s_FBA' % model.getName())
    for (r_id, val) in r_id2val.iteritems():
        r = model.getReaction(r_id)
        if not r:
            r = model.getReaction(format_r_id(r_id, False))
        if not r:
            r = model.getReaction('R_' + format_r_id(r_id, False))
        r.setName('%s: %g' % (r.getId(), val))
    sbml = os.path.join(res_dir, 'Model_FBA.xml')
    libsbml.SBMLWriter().writeSBMLToFile(doc, sbml)
    return sbml


def analyse_by_fba(cobra_model, bm_r_id, directory, objective_sense='maximize', threshold=0,
                   fig_path=None, var_r_id=None, r_ids=None, sbml=None):
    cobra_bm_r_id = format_r_id(bm_r_id)
    optimise_biomass(cobra_model, cobra_bm_r_id, objective_sense, level=logging.DEBUG)
    r_id2val = get_fluxes_larger_than_threshold(cobra_model, threshold=threshold)
    title = '%s reaction %s (%s): %.4g\n==================================\n'\
            % (objective_sense, bm_r_id, get_cobra_r_formula(cobra_model.reactions.get_by_id(cobra_bm_r_id)),
               cobra_model.solution.x_dict[cobra_bm_r_id])
    serialize_fluxes(cobra_model, r_id2val, path=os.path.join(directory, 'fba.txt'), r_ids=r_ids, title=title)

    if fig_path:
        var_r_id = format_r_id(var_r_id)
        var_r = cobra_model.reactions.get_by_id(var_r_id)
        show_reaction_values(cobra_model, var_r_id, {bm_r_id}, bm_r_id,
                             fig_path, "%s flux (micromol/min/gDW) bounds" % var_r.name,
                             "%s production (micromol/min/gDW)" % bm_r_id,
                             "Effect of varying %s on %s" % (var_r.name, bm_r_id),
                             var_r.upper_bound * 2.5, minimized=False, constrained=False)

    if sbml:
        sbml = create_fba_model(sbml, r_id2val, directory)
    return sbml