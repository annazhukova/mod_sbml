from functools import reduce
from itertools import chain
import re

import libsbml

from mod_sbml.sbml.sbml_manager import get_formulas
from mod_sbml.annotation.rdf_annotation_helper import get_is_annotations, get_is_vo_annotations, add_annotation

__author__ = 'anna'

HAS_ROLE_RELATIONSHIP = 'has_role'

COFACTOR_CHEBI_ID = 'chebi:23357'

CONJUGATE_ACID_BASE_RELATIONSHIPS = {'is_conjugate_base_of', 'is_conjugate_acid_of'}

EQUIVALENT_RELATIONSHIPS = CONJUGATE_ACID_BASE_RELATIONSHIPS | {'is_tautomer_of'}

MOLECULAR_ENTITY = 'chebi:23367'

CHEBI_PREFIX = "obo.chebi"

CHEBI_ID_PATTERN = "[cC][Hh][Ee][Bb][Ii]\:\d+"


def get_chebi_id(m):
    for annotation in chain(get_is_annotations(m), get_is_vo_annotations(m)):
        for chebi_id in re.findall(CHEBI_ID_PATTERN, annotation):
            return chebi_id.lower()
    return None


def get_chebi_term_by_annotation(entity, chebi):
    for annotation in chain(get_is_annotations(entity), get_is_vo_annotations(entity)):
        term = chebi.get_term(annotation, check_only_ids=False)
        if term:
            return term
    return None


def infer_chebi_term(m, chebi, model=None):
    term = get_chebi_term_by_annotation(m, chebi)
    if term:
        return term
    names = []
    if model:
        s_type_id = m.getSpeciesType()
        if s_type_id:
            s_type = model.getSpeciesType(s_type_id)
            if s_type:
                term = get_chebi_term_by_annotation(s_type, chebi)
                if term:
                    return term
                names.append(s_type.getName())

    for formula in get_formulas(m):
        if formula and formula != '.':
            term = chebi.get_term(formula, check_only_ids=False)
            if term:
                return term
    name = m.getName()
    names.append(name)
    if name and model:
        c = model.getCompartment(m.getCompartment())
        if c:
            c_name = c.getName() if c.getName() else c.getId()
            name = name.replace('[%s]' % c_name, '').replace(c_name, '').strip()
            names.append(name)
    for name in names:
        if name:
            term = chebi.get_term(name, check_only_ids=False)
            if term:
                return term
    return None


def annotate_metabolites(model, chebi):
    """
    Infers ChEBI terms for metabolites that lack them and annotates.
    :param model: libsbml.Model model of interest
    :param chebi: mod_sbml.onto.obo_ontology.Ontology ChEBI ontology
    :return: void, input model is modified inplace
    """
    for m in model.getListOfSpecies():
        if get_chebi_id(m):
            continue
        term = infer_chebi_term(m, chebi, model)
        if term:
            add_annotation(m, libsbml.BQB_IS, term.get_id(), CHEBI_PREFIX)


def get_species_id2chebi_id(model):
    """
    Gets species id to CheBI term id mapping from the model.
    :param model: libsbml.Model model
    :return: dict {species_id: ChEBI_term_id}
    """
    s_id2chebi_id = {}
    for s in model.getListOfSpecies():
        chebi_id = get_chebi_id(s)
        if chebi_id:
            s_id2chebi_id[s.getId()] = chebi_id
    return s_id2chebi_id


def get_cofactor_ids(onto):
    cofactor_ids = set()
    sub_cofactors = onto.get_descendants(COFACTOR_CHEBI_ID, False)

    def is_cofactor(t_id):
        if COFACTOR_CHEBI_ID == t_id:
            return True
        return onto.get_term(t_id) in sub_cofactors

    for it in onto.get_relationship_participants(HAS_ROLE_RELATIONSHIP):
        subj, rel, obj = it
        if rel == HAS_ROLE_RELATIONSHIP and is_cofactor(obj):
            subj_term = onto.get_term(subj)
            cofactor_ids |= {t.get_id() for t in
                             onto.get_generalized_descendants(subj_term, False, set(), EQUIVALENT_RELATIONSHIPS)}
            cofactor_ids.add(subj_term.get_id())
    return add_equivalent_chebi_ids(onto, cofactor_ids)


def add_equivalent_chebi_ids(ontology, chebi_ids, eq_rel=EQUIVALENT_RELATIONSHIPS):
    """
    Returns a extended set that contains the input ChEBI term id collection
    plus the ids of terms that are equivalent to those in the input collection.
    :param eq_rel: set of equivalence relationships (by default 'is_conjugate_base/acid_of' and 'is_taumer_of').
    :param ontology: mod_sbml.onto.obo_ontology.Ontology ChEBI ontology
    :param chebi_ids: collection of ChEBI term ids.
    :return: set of ChEBI term ids (input + equivalent)
    """
    cup = lambda s1, s2: s1 | s2
    return \
        reduce(cup,
               (reduce(cup,
                       (it.get_all_ids() for it in ontology.get_equivalents(t, relationships=eq_rel)),
                       t.get_all_ids())
                for t in (it for it in (ontology.get_term(ub_id) for ub_id in chebi_ids) if it)), chebi_ids)




