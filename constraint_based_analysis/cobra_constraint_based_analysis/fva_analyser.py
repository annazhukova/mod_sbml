import logging
import os

from constraint_based_analysis.cobra_constraint_based_analysis.fba_manager import serialize_fva, optimise_biomass, \
    get_r_id2fva_bounds
from constraint_based_analysis.cobra_constraint_based_analysis.model_manager import format_r_id
from constraint_based_analysis.efm.efm_serialization_manager import r_ids2sbml
from gibbs.reaction_boundary_manager import set_bounds
from sbml.sbml_manager import reverse_reaction
from serialization.serialization_manager import get_cobra_r_formula

__author__ = 'anna'



def analyse_by_fva(cobra_model, bm_r_id, directory, objective_sense='maximize', threshold=0, r_ids=None, sbml=None,
                   rewrite=True):
    fva_file = os.path.join(directory, 'fva.txt')
    if not rewrite and os.path.exists(fva_file):
        if sbml:
            new_sbml = os.path.join(directory, 'Model_FVA.xml')
            if os.path.exists(new_sbml):
                return new_sbml
        else:
            return sbml
    cobra_bm_r_id = format_r_id(bm_r_id)
    optimise_biomass(cobra_model, cobra_bm_r_id, objective_sense, level=logging.DEBUG)
    r_id2bounds = get_r_id2fva_bounds(cobra_model, threshold=threshold)
    title = '%s reaction %s (%s): %.4g\n==================================\n'\
            % (objective_sense, bm_r_id, get_cobra_r_formula(cobra_model.reactions.get_by_id(cobra_bm_r_id)),
               cobra_model.solution.x_dict[cobra_bm_r_id])
    serialize_fva(cobra_model, r_id2bounds, fva_file, r_ids=r_ids, title=title)

    if sbml:
        new_sbml = os.path.join(directory, 'Model_FVA.xml')
        essential_r_id2rev = {r_id: u < 0 for (r_id, (l, u)) in r_id2bounds.iteritems() if l * u > 0}
        create_essential_r_ids_model(sbml, essential_r_id2rev, directory)
        sbml = create_fva_model(sbml, r_id2bounds, new_sbml)
    return r_id2bounds, sbml


def create_fva_model(sbml, r_id2bounds, new_sbml):
    def r_updater(r):
        l_b, u_b = r_id2bounds[format_r_id(r.id)]
        if u_b < 0:
            reverse_reaction(r)
            r.setName('-%s' % r.id)
            set_bounds(r, -u_b, -l_b)
        else:
            r.setName(r.getId())
            set_bounds(r, l_b, u_b)
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

    r_ids = set(r_id2rev.iterkeys()) | {format_r_id(r_id, False) for r_id in r_id2rev.iterkeys()}
    new_sbml = os.path.join(res_dir, 'Model_common.xml')
    r_ids2sbml(r_ids, sbml, new_sbml, 'common', r_updater)
    return new_sbml