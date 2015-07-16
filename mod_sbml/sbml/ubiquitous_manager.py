from itertools import chain

from mod_sbml.annotation.chebi.chebi_annotator import get_species_to_chebi, get_cofactor_ids, \
    add_equivalent_chebi_ids
from mod_sbml.onto.onto_getter import get_chebi
from mod_sbml.onto.obo_ontology import parse
from mod_sbml.sbml.sbml_manager import get_reactants, get_products

__author__ = 'anna'

UBIQUITOUS_THRESHOLD = 20


# most common ones, like water, H+, oxygen, NAD, etc.
COMMON_UB_IDS = {'chebi:17808', 'chebi:37568', 'chebi:24636', 'chebi:15422', 'chebi:15377', 'chebi:15378',
                 'chebi:15379', 'chebi:37563', 'chebi:37565', 'chebi:17877', 'chebi:18245', 'chebi:16311',
                 'chebi:16192', 'chebi:15846', 'chebi:30089', 'chebi:16497', 'chebi:25805', 'chebi:26020',
                 'chebi:13390', 'chebi:13392', 'chebi:43474', 'chebi:57783', 'chebi:17625', 'chebi:18359',
                 'chebi:25524', 'chebi:58115', 'chebi:16695', 'chebi:15346', 'chebi:17552', 'chebi:58107',
                 'chebi:17330', 'chebi:57692', 'chebi:17013', 'chebi:61404', 'chebi:30616', 'chebi:61402',
                 'chebi:57299', 'chebi:58307', 'chebi:58223', 'chebi:35780', 'chebi:17659', 'chebi:28862',
                 'chebi:58342', 'chebi:57287', 'chebi:17093', 'chebi:16027', 'chebi:58189', 'chebi:33813',
                 'chebi:29888', 'chebi:16908', 'chebi:57540', 'chebi:456216', 'chebi:456215', 'chebi:57945',
                 'chebi:28971', 'chebi:16474', 'chebi:28850', 'chebi:18421', 'chebi:16761', 'chebi:17361',
                 'chebi:26689', 'chebi:57288', 'chebi:16284', 'chebi:16234', 'chebi:16174', 'chebi:82680',
                 'chebi:18075', 'chebi:16171', 'chebi:18077', 'chebi:15713', 'chebi:16039', 'chebi:16238',
                 'chebi:36080', 'chebi:29375', 'chebi:73342', 'chebi:15366', 'chebi:16526', 'chebi:17984',
                 'chebi:17544', 'chebi:28846', 'chebi:29033', 'chebi:17677', 'chebi:15996', 'chebi:18009',
                 'chebi:58210', 'chebi:61429', 'chebi:58280', 'chebi:18361', 'chebi:15918', 'chebi:246422',
                 'chebi:17239', 'chebi:73627', 'chebi:13389', 'chebi:16240', 'chebi:58245', 'chebi:29325',
                 'chebi:13534', 'chebi:15351', 'chebi:26078'}


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
                             species_id2chebi_id=None, model=None, onto=None):
    if not onto:
        onto = parse(get_chebi())
    ubiquitous_chebi_ids = set()
    if add_cofactors:
        ubiquitous_chebi_ids |= get_cofactor_ids(onto)
    if add_common:
        ubiquitous_chebi_ids = add_equivalent_chebi_ids(onto, COMMON_UB_IDS)
    if add_frequent and model:
        if not species_id2chebi_id:
            species_id2chebi_id = get_species_to_chebi(model, onto)
        ubiquitous_chebi_ids |= \
            add_equivalent_chebi_ids(onto, get_frequent_term_ids(model, species_id2chebi_id, ubiquitous_threshold))
    return ubiquitous_chebi_ids


def select_metabolite_ids_by_term_ids(model, selected_chebi_ids, species_id2chebi_id=None):
    if not species_id2chebi_id:
        onto = parse(get_chebi())
        species_id2chebi_id = get_species_to_chebi(model, onto)
    return {s.getId() for s in model.getListOfSpecies() if
            s.getId() in species_id2chebi_id and species_id2chebi_id[s.getId()] in selected_chebi_ids}
