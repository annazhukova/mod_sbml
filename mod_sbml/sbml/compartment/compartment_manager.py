from collections import defaultdict
from itertools import chain
import logging
from mod_sbml.annotation.chebi.chebi_annotator import get_chebi_id

from mod_sbml.sbml.sbml_manager import get_products, create_species, get_reactants, create_reaction, \
    create_compartment

BOUNDARY_C_NAME = 'Boundary'
BOUNDARY_C_ID = 'Boundary'

__author__ = 'anna'


def separate_boundary_species(model):
    """
    Creates a boundary compartment with the id 'Boundary' and moves the boundary species there.
    :param model: object of libsbml.Model
    :return: void
    """
    boundary_comp = create_boundary_compartment_if_needed(model)

    key2boundary_s_id = {}
    for r in model.getListOfReactions():
        rs = {it for it in (model.getSpecies(species_ref.getSpecies()) for species_ref in r.getListOfReactants()) if it}
        ps = {it for it in (model.getSpecies(species_ref.getSpecies()) for species_ref in r.getListOfProducts()) if it}
        chebi_id2ss = defaultdict(list)
        for s in chain(rs, ps):
            chebi_id = get_chebi_id(s)
            if chebi_id:
                chebi_id2ss[chebi_id].append(s)
        for s in (s for s in chain(rs, ps) if s.getBoundaryCondition()):
            s_id = s.getId()
            chebi_id = get_chebi_id(s)
            if chebi_id and len(chebi_id2ss[chebi_id]) > 1 \
                    and len({it.getCompartment() for it in chebi_id2ss[chebi_id]}) == 1:
                s.setCompartment(boundary_comp.getId())
                key2boundary_s_id[chebi_id] = s_id
            else:
                if not chebi_id:
                    chebi_id = s_id
                if chebi_id in key2boundary_s_id:
                    boundary_s_id = key2boundary_s_id[chebi_id]
                    if boundary_s_id == s_id:
                        continue
                else:
                    boundary_s_id = \
                        create_species(model, compartment_id=boundary_comp.getId(), name=s.getName(), bound=True,
                                       id_='%s_b' % s.getId(), type_id=s.getSpeciesType(),
                                       sbo_id=s.getSBOTerm()).getId()
                    key2boundary_s_id[chebi_id] = boundary_s_id
                create_reaction(model, {boundary_s_id: 1}, {s_id: 1},
                                name='Exchange %s' % (s.getName() if s.getName() else s_id),
                                reversible=True, id_='%s_exchange' % s_id)
                s.setBoundaryCondition(False)
    create_boundary_species_in_boundary_reactions(model, key2boundary_s_id, boundary_comp)


def create_boundary_compartment_if_needed(model):
    boundary_comp = get_boundary_compartment(model)
    if boundary_comp:
        return boundary_comp
    return create_compartment(model, name=BOUNDARY_C_NAME, id_=BOUNDARY_C_ID)


def create_boundary_species_in_boundary_reactions(model, key2boundary_s_id=None, boundary_comp=None):
    if not boundary_comp:
        boundary_comp = create_boundary_compartment_if_needed(model)
    if not key2boundary_s_id:
        key2boundary_s_id = {}
    for r in model.getListOfReactions():
        if r.getNumReactants() == 0:
            for s_id, st in get_products(r, stoichiometry=True):
                species = model.getSpecies(s_id)
                if not species:
                    logging.error('Check your model: reaction %s has an undefined product %s' % (r.getId(), s_id))
                    continue
                key = get_chebi_id(species)
                if not key:
                    key = s_id
                if key in key2boundary_s_id:
                    boundary_s_id = key2boundary_s_id[key]
                else:
                    boundary_s_id = \
                        create_species(model, compartment_id=boundary_comp.getId(), name=species.getName(), bound=True,
                                       id_='%s_b' % species.getId(), type_id=species.getSpeciesType(),
                                       sbo_id=species.getSBOTerm()).getId()
                    key2boundary_s_id[key] = boundary_s_id
                new_m = r.createReactant()
                new_m.setSpecies(boundary_s_id)
                new_m.setStoichiometry(st)
        if r.getNumProducts() == 0:
            for s_id, st in get_reactants(r, stoichiometry=True):
                species = model.getSpecies(s_id)
                if not species:
                    logging.error('Check your model: reaction %s has an undefined reactant %s' % (r.getId(), s_id))
                    continue
                key = get_chebi_id(species)
                if not key:
                    key = s_id
                if key in key2boundary_s_id:
                    boundary_s_id = key2boundary_s_id[key]
                else:
                    boundary_s_id = \
                        create_species(model, compartment_id=boundary_comp.getId(), name=species.getName(), bound=True,
                                       id_='%s_b' % species.getId(), type_id=species.getSpeciesType(),
                                       sbo_id=species.getSBOTerm()).getId()
                    key2boundary_s_id[key] = boundary_s_id
                new_m = r.createProduct()
                new_m.setSpecies(boundary_s_id)
                new_m.setStoichiometry(st)


def get_boundary_compartment(model):
    for compartment in model.getListOfCompartments():
        if compartment.getName() and BOUNDARY_C_NAME.lower() == compartment.getName().lower() \
                or BOUNDARY_C_ID.lower() == compartment.getId().lower() or 'b' == compartment.getId().lower():
            return compartment
    return None


def need_boundary_compartment(model):
    """
    Checks if the model does not contain a Boundary compartment yet but contains at least one reaction
    using metabolites with the same CHEBI id and the same compartment, but different boundary condition (True and False)
    :param model: object of libsbml.Model
    :return: if the model would benefit from the creation of a boundary compartment
    """
    b_comp = get_boundary_compartment(model)
    if b_comp:
        return False
    for r in model.getListOfReactions():
        if r.getNumReactants() == 0 or r.getNumProducts() == 0:
            return True
        if (r.getNumReactants() > 1 or r.getNumProducts() > 1) \
                and next((it for it in chain(get_reactants(r), get_products(r))
                          if model.getSpecies(it).getBoundaryCondition()), False):
            return True
    return False
