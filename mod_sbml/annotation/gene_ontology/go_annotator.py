from itertools import chain
import re

import libsbml

from mod_sbml.annotation.rdf_annotation_helper import add_annotation, get_is_annotations, get_is_vo_annotations

GO_PREFIX = 'go'

__author__ = 'anna'

GO_ID_PATTERN = "[gG][Oo]\:\d+"


def get_go_id(c):
    for annotation in chain(get_is_annotations(c), get_is_vo_annotations(c)):
        for go_id in re.findall(GO_ID_PATTERN, annotation):
            return go_id.lower()
    return None


def infer_go_term(comp, onto):
    """
    Find a Gene Ontology term for a compartment of interest,
    using its name.
    :param comp: libSBML Compartment
    :param onto: the Gene Ontology
    :return: term if it was found otherwise None
    """
    for annotation in chain(get_is_annotations(comp), get_is_vo_annotations(comp)):
        term = onto.get_term(annotation, check_only_ids=False)
        if term:
            return term
    return onto.get_term(comp.getName() if comp.getName() else comp.getId(), check_only_ids=False)


def annotate_compartments(model, go):
    for comp in model.getListOfCompartments():
        if get_go_id(comp):
            continue
        term = infer_go_term(comp, go)
        if term:
            add_annotation(comp, libsbml.BQB_IS, term.get_id(), GO_PREFIX)
