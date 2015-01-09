from collections import defaultdict
import os
from libsbml import Species, BQB_IS, BQB_IS_VERSION_OF
import misc
from obo_ontology import miriam_to_term_id, to_identifiers_org_format
from sbml_manager import get_qualifier_values, add_annotation

__author__ = 'anna'


def get_chebi():
    return "%s/data/chebi.obo" % os.path.dirname(os.path.abspath(misc.__file__))


def get_is_annotations(entity):
    return (miriam_to_term_id(it) for it in get_qualifier_values(entity.getAnnotation(), BQB_IS))


def get_is_vo_annotations(entity):
    return (miriam_to_term_id(it) for it in get_qualifier_values(entity.getAnnotation(), BQB_IS_VERSION_OF))


def get_term(entity, chebi):
    for is_annotation in get_is_annotations(entity):
        term = chebi.get_term(is_annotation)
        if term:
            return term
    for is_vo_annotation in get_is_vo_annotations(entity):
        term = chebi.get_term(is_vo_annotation)
        if term:
            return term
    return None


def normalize(name):
    return name.strip()


def get_names(entity):
    name = normalize(entity.getName())
    name_bis = name
    end = name_bis.find("(")
    if end != -1 and end != 0:
        name_bis = name_bis[0:end].strip()
    return name, name_bis


def get_species_term(species, chebi, model):
    term = get_term(species, chebi)
    if not term:
        s_type_id = species.getSpeciesType()
        if s_type_id:
            s_type = model.getSpeciesType(s_type_id)
            if s_type:
                term = get_term(s_type, chebi)
    return term


def get_species_to_chebi(model, chebi, guess=True):
    species2chebi = {}
    used_terms = set()
    entity2species = defaultdict(set)

    # process annotated ones
    # and find those that need to be annotated
    for species in model.getListOfSpecies():
        term = get_term(species, chebi)
        entity = species
        if not term:
            s_type_id = species.getSpeciesType()
            if s_type_id:
                s_type = model.getSpeciesType(s_type_id)
                if s_type:
                    entity = s_type
                    term = get_term(s_type, chebi)
        if term:
            species2chebi[species.getId()] = term.get_id()
            used_terms.add(term)
            continue
        else:
            entity2species[entity].add(species)
    if guess:
        # annotate unannotated
        for entity, species_set in entity2species.iteritems():
            name, name_bis = get_names(entity)
            if isinstance(entity, Species):
                index = name.lower().find(
                    "[{0}]".format(model.getCompartment(entity.getCompartment()).getName().lower()))
                if -1 != index:
                    name = name[:index].strip()
            possibilities = chebi.get_ids_by_name(name)
            if not possibilities:
                possibilities = chebi.get_ids_by_name(name_bis)
            if not possibilities:
                continue
            possibilities = {chebi.get_term(it) for it in possibilities}
            intersection = possibilities & used_terms
            term = intersection.pop() if intersection else possibilities.pop()
            for species in species_set:
                species2chebi[species.getId()] = term.get_id()
            add_annotation(entity, BQB_IS, to_identifiers_org_format(term.get_id(), "obo.chebi"))
            used_terms.add(term)
    return species2chebi
