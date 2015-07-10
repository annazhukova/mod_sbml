import logging

from mod_sbml.constraint_based_analysis.cobra_constraint_based_analysis.fba_manager import optimise_biomass, \
    get_fluxes_larger_than_threshold, round_value
from mod_sbml.constraint_based_analysis.efm.efm_serialization_manager import r_ids2sbml
from mod_sbml.constraint_based_analysis.cobra_constraint_based_analysis.model_manager import format_r_id
from mod_sbml.sbml.sbml_manager import reverse_reaction

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


def analyse_by_fba(cobra_model, bm_r_id, objective_sense='maximize', threshold=0):
    cobra_bm_r_id = format_r_id(bm_r_id)
    opt_value = optimise_biomass(cobra_model, cobra_bm_r_id, objective_sense, level=logging.DEBUG)
    opt_value = round_value(opt_value)
    r_id2val = get_fluxes_larger_than_threshold(cobra_model, threshold=threshold)

    return r_id2val, opt_value