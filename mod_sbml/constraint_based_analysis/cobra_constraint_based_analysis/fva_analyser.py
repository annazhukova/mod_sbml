import logging
import os

from mod_sbml.constraint_based_analysis.cobra_constraint_based_analysis.fba_manager import optimise_biomass, \
    get_r_id2fva_bounds, round_value
from mod_sbml.constraint_based_analysis.efm.efm_serialization_manager import r_ids2sbml
from mod_sbml.constraint_based_analysis.cobra_constraint_based_analysis.model_manager import format_r_id
from mod_sbml.gibbs.reaction_boundary_manager import set_bounds
from mod_sbml.sbml.sbml_manager import reverse_reaction

__author__ = 'anna'



def analyse_by_fva(cobra_model, bm_r_id, objective_sense='maximize', threshold=0):
    cobra_bm_r_id = format_r_id(bm_r_id)
    opt_value = optimise_biomass(cobra_model, cobra_bm_r_id, objective_sense, level=logging.DEBUG)
    opt_value = round_value(opt_value)
    r_id2bounds = get_r_id2fva_bounds(cobra_model, threshold=threshold)

    return r_id2bounds, opt_value


def create_fva_model(sbml, r_id2bounds, new_sbml):
    def r_updater(r):
        l_b, u_b = r_id2bounds[format_r_id(r.id)]
        # if u_b < 0:
        #     reverse_reaction(r)
        #     r.setName('-%s' % r.id)
        #     set_bounds(r, -u_b, -l_b)
        # else:
        #     r.setName(r.getId())
        #     set_bounds(r, l_b, u_b)
        set_bounds(r, l_b, max(u_b, 0))
        if l_b * u_b >= 0:
            r.setReversible(False)

    r_ids = set(r_id2bounds.iterkeys()) | {format_r_id(r_id, False) for r_id in r_id2bounds.iterkeys()}
    r_ids2sbml(r_ids, sbml, new_sbml, 'FVA', r_updater)
    return new_sbml


def create_essential_r_ids_model(sbml, r_id2rev, res_dir):
    def r_updater(r):
        rev = r_id2rev[format_r_id(r.id)]
        if rev < 0:
            reverse_reaction(r)
            r.setName('-%s' % r.id)
        else:
            r.setName(r.getId())
        r.setReversible(False)

    r_ids = set(r_id2rev.iterkeys()) | {format_r_id(r_id, False) for r_id in r_id2rev.iterkeys()}
    new_sbml = os.path.join(res_dir, 'Model_essential_reactions.xml')
    r_ids2sbml(r_ids, sbml, new_sbml, 'essential_reactions', r_updater)
    return new_sbml