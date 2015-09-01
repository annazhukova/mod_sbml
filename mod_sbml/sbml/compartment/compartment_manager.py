from collections import defaultdict
from itertools import chain
import logging

import libsbml

from mod_sbml.sbml.sbml_manager import generate_unique_id

BOUNDARY_C_NAME = 'Boundary'

BOUNDARY_C_ID = 'Boundary'

__author__ = 'anna'


def separate_boundary_species(model):
    """
    Creates a boundary compartment with the id 'Boundary' and moves the boundary species there.
    :param model: object of libsbml.Model
    :return: void
    """
    boundary_comp = model.createCompartment()
    id_ = generate_unique_id(model, BOUNDARY_C_ID)
    if libsbml.LIBSBML_OPERATION_SUCCESS != boundary_comp.setId(id_):
        logging.error("boundary compartment  %s creation error" % id_)
    boundary_comp.setName(BOUNDARY_C_NAME)
    for s in model.getListOfSpecies():
        if s.getBoundaryCondition():
            s.setCompartment(boundary_comp.getId())


def need_boundary_compartment(model, s_id2chebi):
    """
    Checks if the model does not contain a Boundary compartment yet but contains at least one reaction
    using metabolites with the same CHEBI id and the same compartment, but different boundary condition (True and False).
    :param model: object of libsbml.Model
    :param s_id2chebi: mapping between metabolite ids in the model and their CHEBI ids
    :return: if the model would benefit from the creation of a boundary compartment
    """
    for compartment in model.getListOfCompartments():
        if compartment.getName() and BOUNDARY_C_NAME.lower() == compartment.getName().lower():
            return False
    for r in model.getListOfReactions():
        elements = chain((model.getSpecies(species_ref.getSpecies()) for species_ref in r.getListOfReactants()),
                         (model.getSpecies(species_ref.getSpecies()) for species_ref in r.getListOfProducts()))
        chebi_id2ss = defaultdict(list)
        for s in elements:
            if not s:
                raise ValueError(
                    "Your model includes undeclared metabolites in the reaction %s. Please, fix your model." % r.getid())
            if s.getId() in s_id2chebi:
                chebi_id2ss[s_id2chebi[s.getId()]].append(s)
        for ss in chebi_id2ss.itervalues():
            if len(ss) > 1:
                for s_b in (s for s in ss if s.getBoundaryCondition()):
                    if next((s for s in ss
                             if not s.getBoundaryCondition() and s_b.getCompartment() == s.getCompartment()), None):
                        return True
    return False
