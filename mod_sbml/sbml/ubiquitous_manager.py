from itertools import chain

from mod_sbml.annotation.chebi.chebi_annotator import get_species_to_chebi, add_equivalent_chebi_ids
from mod_sbml.annotation.chebi.chebi_serializer import get_chebi, COMMON_TERMS_FILE, COFACTORS_FILE
from mod_sbml.onto import parse_simple
from mod_sbml.sbml.sbml_manager import get_reactants, get_products

__author__ = 'anna'

UBIQUITOUS_THRESHOLD = 20

common_ch_ids = None
cofactor_ch_ids = None


def get_frequent_term_ids(model, species_id2chebi_id, threshold=UBIQUITOUS_THRESHOLD):
    """
    The function returns a set of identifiers of ChEBI term ids of ubiquitous metabolites belonging to the given model.
    The species in the model are divided into two groups: ubiquitous ones and the others.
    Ubiquitous species are those participating in more than {@link #threshold threshold number} of reactions.

    :param model: a {@link #libsbml.Model Model} object.
    :param species_id2chebi_id: a mapping between species identifiers (string) and their ChEBI identifiers (string).
    :param threshold: (Optional) A minimal number of reactions a species should participate in to become a ubiquitous one.
    The default value is {@link #UBIQUITOUS_THRESHOLD UBIQUITOUS_THRESHOLD}.
    :return: A set of ubiquitous ChEBI term identifiers.
    """
    chebi2vote = {}

    for reaction in model.getListOfReactions():
        for element in chain(get_reactants(reaction), get_products(reaction)):
            # if we do not have a ChEBI ubiquitous for it,
            # it will be considered ubiquitous anyway
            if element not in species_id2chebi_id:
                continue
            chebi_id = species_id2chebi_id[element]
            compartment = model.getSpecies(element).getCompartment()
            key = chebi_id, compartment
            if key in chebi2vote:
                chebi2vote[key] += 1
            else:
                chebi2vote[key] = 1

    ubiquitous_ch_ids = set()
    for element, vote in chebi2vote.iteritems():
        if vote > threshold:
            ubiquitous_ch_ids.add(element[0])

    return ubiquitous_ch_ids


def get_ubiquitous_chebi_ids(add_common=True, add_cofactors=True, add_frequent=False,
                             ubiquitous_threshold=UBIQUITOUS_THRESHOLD,
                             species_id2chebi_id=None, model=None, chebi=None):
    global cofactor_ch_ids, common_ch_ids
    ubiquitous_chebi_ids = set()
    if add_cofactors:
        if not cofactor_ch_ids:
            with open(COFACTORS_FILE, 'r') as f:
                cofactor_ch_ids = set(f.readline().split('\t'))
        ubiquitous_chebi_ids |= cofactor_ch_ids
    if add_common:
        if not common_ch_ids:
            with open(COMMON_TERMS_FILE, 'r') as f:
                common_ch_ids = set(f.readline().split('\t'))
        ubiquitous_chebi_ids = common_ch_ids
    if add_frequent and model:
        if not chebi:
            chebi = parse_simple(get_chebi())
        if not species_id2chebi_id:
            species_id2chebi_id = get_species_to_chebi(model, chebi)
        ubiquitous_chebi_ids |= \
            add_equivalent_chebi_ids(chebi, get_frequent_term_ids(model, species_id2chebi_id, ubiquitous_threshold))
    return ubiquitous_chebi_ids


def select_metabolite_ids_by_term_ids(model, selected_chebi_ids, species_id2chebi_id=None):
    if not species_id2chebi_id:
        chebi = parse_simple(get_chebi())
        species_id2chebi_id = get_species_to_chebi(model, chebi)
    return {s.getId() for s in model.getListOfSpecies() if
            s.getId() in species_id2chebi_id and species_id2chebi_id[s.getId()] in selected_chebi_ids}

