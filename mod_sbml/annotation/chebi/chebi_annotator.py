from collections import defaultdict

import libsbml

from mod_sbml.annotation.miriam_converter import to_identifiers_org_format
from mod_sbml.sbml.sbml_manager import get_formulas
from mod_sbml.annotation.rdf_annotation_helper import get_is_annotations, get_is_vo_annotations, get_annotations, \
    add_annotation

__author__ = 'anna'

HAS_ROLE_RELATIONSHIP = 'has_role'

COFACTOR_CHEBI_ID = 'chebi:23357'

CONJUGATE_ACID_BASE_RELATIONSHIPS = {'is_conjugate_base_of', 'is_conjugate_acid_of'}

EQUIVALENT_RELATIONSHIPS = CONJUGATE_ACID_BASE_RELATIONSHIPS | {'is_tautomer_of'}

MOLECULAR_ENTITY = 'chebi:23367'

CHEBI_PREFIX = "obo.chebi"


def get_chebi_id(m):
    for annotation in get_annotations(m, libsbml.BQB_IS):
        if annotation.lower().find("chebi") != -1:
            return annotation.lower().replace("obo.chebi:", '').replace("chebi:chebi:", 'chebi:')
    return None


def get_term(entity, chebi):
    for is_annotation in get_is_annotations(entity):
        term = chebi.get_term(is_annotation, check_only_ids=False)
        if term:
            return term
    for is_vo_annotation in get_is_vo_annotations(entity):
        term = chebi.get_term(is_vo_annotation, check_only_ids=False)
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


def find_term_id(entity, chebi):
    term = get_term(entity, chebi)
    if term:
        return term.get_id()

    for formula in get_formulas(entity):
        if formula and formula != '.':
            term = chebi.get_term(formula, check_only_ids=False)
            if term:
                return term.get_id()
    return None


def get_species_to_chebi(model, chebi, guess=True):
    species2chebi = {}
    s_type_id2chebi = {}
    unannotated = []

    # process species types
    for s_type in model.getListOfSpeciesTypes():
        t_id = find_term_id(s_type, chebi)
        if t_id:
            if not next((annotation for annotation in get_is_annotations(s_type) if annotation.find('chebi') != -1),
                         False):
                add_annotation(s_type, libsbml.BQB_IS, to_identifiers_org_format(t_id, CHEBI_PREFIX))
            s_type_id2chebi[s_type.getId()] = t_id

    # process species
    for species in model.getListOfSpecies():
        s_type_id = species.getSpeciesType()
        if s_type_id and s_type_id in s_type_id2chebi:
            t_id = s_type_id2chebi[s_type_id]
        else:
            t_id = find_term_id(species, chebi)
        if t_id:
            species2chebi[species.getId()] = t_id
            if not next((annotation for annotation in get_is_annotations(species) if annotation.find('chebi') != -1),
                        False):
                add_annotation(species, libsbml.BQB_IS, to_identifiers_org_format(t_id, CHEBI_PREFIX))
            if s_type_id and s_type_id not in s_type_id2chebi:
                s_type_id2chebi[s_type_id] = t_id
                add_annotation(model.getSpeciesType(s_type_id), libsbml.BQB_IS,
                               to_identifiers_org_format(t_id, CHEBI_PREFIX))
        else:
            unannotated.append(species)

    s_t_id2unannotated = defaultdict(list)
    for species in unannotated:
        s_type_id = species.getSpeciesType()
        if s_type_id and s_type_id in s_type_id2chebi:
            t_id = s_type_id2chebi[s_type_id]
            species2chebi[species.getId()] = t_id
            add_annotation(species, libsbml.BQB_IS, to_identifiers_org_format(t_id, CHEBI_PREFIX))
        else:
            if s_type_id:
                s_t_id2unannotated[s_type_id].append(species)
            else:
                s_t_id2unannotated[species.getId()].append(species)

    # annotate unannotated
    if guess:
        for s_t_id, species_list in s_t_id2unannotated.iteritems():
            s_type = model.getSpeciesType(s_t_id)
            name = ''
            if s_type:
                name = normalize(s_type.getName())
            for species in species_list:
                if name:
                    break
                c_name = normalize("[{0}]".format(model.getCompartment(species.getCompartment()).getName()))
                name = normalize(species.getName()).replace(c_name, '').strip()
            if not name:
                continue
            term = chebi.get_term(name, check_only_ids=False)
            if not term:
                continue
            t_id = term.get_id()
            for species in species_list:
                species2chebi[species.getId()] = t_id
                add_annotation(species, libsbml.BQB_IS, to_identifiers_org_format(t_id, CHEBI_PREFIX))
            if s_type:
                add_annotation(s_type, libsbml.BQB_IS, to_identifiers_org_format(t_id, CHEBI_PREFIX))
    return species2chebi


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


def add_equivalent_chebi_ids(ontology, chebi_ids):
    return reduce(lambda s1, s2: s1 | s2,
                  (reduce(lambda s1, s2: s1 | s2,
                          (it.get_all_ids() for it in ontology.get_equivalents(t, relationships=EQUIVALENT_RELATIONSHIPS)),
                          t.get_all_ids())
                   for t in (it for it in (ontology.get_term(ub_id) for ub_id in chebi_ids) if it)), chebi_ids)




