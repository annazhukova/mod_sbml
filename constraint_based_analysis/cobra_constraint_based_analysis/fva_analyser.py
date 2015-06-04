import logging
import os
import libsbml
from constraint_based_analysis.cobra_constraint_based_analysis.fba_manager import serialize_fva, optimise_biomass, \
    get_r_id2fva_bounds
from constraint_based_analysis.cobra_constraint_based_analysis.model_manager import format_r_id
from gibbs.reaction_boundary_manager import set_bounds
from sbml.sbml_manager import submodel, reverse_reaction
from serialization.serialization_manager import get_cobra_r_formula

__author__ = 'anna'



def analyse_by_fva(cobra_model, bm_r_id, directory, objective_sense='maximize', threshold=0, r_ids=None, sbml=None):
    cobra_bm_r_id = format_r_id(bm_r_id)
    optimise_biomass(cobra_model, cobra_bm_r_id, objective_sense, level=logging.DEBUG)
    r_id2bounds = get_r_id2fva_bounds(cobra_model, threshold=threshold)
    title = '%s reaction %s (%s): %.4g\n==================================\n'\
            % (objective_sense, bm_r_id, get_cobra_r_formula(cobra_model.reactions.get_by_id(cobra_bm_r_id)),
               cobra_model.solution.x_dict[cobra_bm_r_id])
    serialize_fva(cobra_model, r_id2bounds, os.path.join(directory, 'fva.txt'), r_ids=r_ids, title=title)

    if sbml:
        essential_r_id2rev = {r_id: u < 0 for (r_id, (l, u)) in r_id2bounds.iteritems() if l * u > 0}
        create_essential_r_ids_model(sbml, essential_r_id2rev, directory)
        sbml = create_fva_model(sbml, r_id2bounds, directory)
    return sbml


def create_fva_model(sbml, r_id2bounds, res_dir):
    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()

    sbml = os.path.join(res_dir, 'Model_FVA.xml')
    submodel(set(r_id2bounds.iterkeys()) | {format_r_id(r_id, False) for r_id in r_id2bounds.iterkeys()} |
             {'R_' + format_r_id(r_id, False) for r_id in r_id2bounds.iterkeys()}, model)
    model.setId('%s_FVA' % model.getId())
    model.setName('%s_FVA' % model.getName())
    for r_id, (l_b, u_b) in r_id2bounds.iteritems():
        r = model.getReaction(r_id)
        if not r:
            r = model.getReaction(format_r_id(r_id, False))
        if not r:
            r = model.getReaction('R_' + format_r_id(r_id, False))
        set_bounds(r, l_b, u_b)
    libsbml.SBMLWriter().writeSBMLToFile(doc, sbml)
    return sbml


def create_essential_r_ids_model(sbml, r_id2rev, res_dir):
    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()
    m_id = model.getId()
    m_name = model.getName()

    sbml = os.path.join(res_dir, 'Model_common.xml')
    submodel(set(r_id2rev.iterkeys()) | {format_r_id(r_id, False) for r_id in r_id2rev.iterkeys()} |
             {'R_' + format_r_id(r_id, False) for r_id in r_id2rev.iterkeys()}, model)
    model.setId('%s_common' % m_id)
    model.setName('%s_common' % m_name)
    for r_id, rev in r_id2rev.iteritems():
        if rev:
            r = model.getReaction(r_id)
            if not r:
                r = model.getReaction(format_r_id(r_id, False))
            if not r:
                r = model.getReaction('R_' + format_r_id(r_id, False))
            reverse_reaction(r)
    libsbml.SBMLWriter().writeSBMLToFile(doc, sbml)
    return sbml