from mod_sbml.annotation.chebi.chebi_annotator import get_species_to_chebi, get_cofactors, \
    add_equivalent_chebi_ids

from mod_sbml.onto.onto_getter import get_chebi
from mod_sbml.onto.obo_ontology import parse
from mod_sbml.sbml.sbml_manager import get_metabolites

__author__ = 'anna'

UBIQUITOUS_THRESHOLD = 20


# most common ones, like water, H+, oxygen, NAD, etc.
COMMON_UB_IDS = {'chebi:37568', 'chebi:57783', 'chebi:17625', 'chebi:37563', 'chebi:17552', 'chebi:17361',
                 'chebi:16311', 'chebi:16192', 'chebi:15846', 'chebi:61429', 'chebi:16234', 'chebi:16174',
                 'chebi:58210', 'chebi:16171', 'chebi:36080', 'chebi:15713', 'chebi:16238', 'chebi:43474',
                 'chebi:15378', 'chebi:15379', 'chebi:58115', 'chebi:29375', 'chebi:16695', 'chebi:58342',
                 'chebi:15346', 'chebi:37565', 'chebi:16526', 'chebi:17544', 'chebi:17013', 'chebi:61404',
                 'chebi:30616', 'chebi:18009', 'chebi:58307', 'chebi:58223', 'chebi:18361', 'chebi:28862',
                 'chebi:15918', 'chebi:246422', 'chebi:28850', 'chebi:16240', 'chebi:58245', 'chebi:16908',
                 'chebi:13534', 'chebi:456216', 'chebi:456215', 'chebi:15351', 'chebi:30089', 'chebi:15422',
                 'chebi:57299', 'chebi:25805', 'chebi:26689', 'chebi:13390', 'chebi:57540', 'chebi:25524',
                 'chebi:13389', 'chebi:13392', 'chebi:28971', 'chebi:17984', 'chebi:29888', 'chebi:26020',
                 'chebi:73342', 'chebi:35780', 'chebi:26078', 'chebi:24636'} \
                | {'chebi:15422', 'chebi:15846', 'chebi:15378', 'chebi:16908', 'chebi:16027', 'chebi:16474',
                   'chebi:16761', 'chebi:17361', 'chebi:15713', 'chebi:17877', 'chebi:15366', 'chebi:17544',
                   'chebi:17677', 'chebi:16240', 'chebi:15377', 'chebi:15379', 'chebi:17552', 'chebi:16311',
                   'chebi:15346', 'chebi:17659', 'chebi:28862', 'chebi:16238', 'chebi:16526', 'chebi:17984',
                   'chebi:16174', 'chebi:18009', 'chebi:15351', 'chebi:16039', 'chebi:18421', 'chebi:29033',
                   'chebi:18075', 'chebi:18077', 'chebi:17808', 'chebi:18359', 'chebi:16497', 'chebi:16284',
                   'chebi:28846', 'chebi:15996', 'chebi:17239', 'chebi:37565', 'chebi:18245', 'chebi:57287',
                   'chebi:73342', 'chebi:33813', 'chebi:57783', 'chebi:57945', 'chebi:29375', 'chebi:82680',
                   'chebi:58280', 'chebi:30616', 'chebi:61402', 'chebi:58189', 'chebi:57288', 'chebi:57692',
                   'chebi:58307', 'chebi:28971', 'chebi:17093', 'chebi:456215', 'chebi:13534', 'chebi:73627',
                   'chebi:17330', 'chebi:58107', 'chebi:29325'}


def get_ubiquitous_species_set(model, species_id2chebi_id, ontology, threshold=UBIQUITOUS_THRESHOLD):
    """
    The function returns a set of identifiers of ubiquitous species belonging to the given model.
    The species in the model are divided into two groups: ubiquitous ones and the others.
    Ubiquitous species are those participating in more than {@link #threshold threshold number} of reactions.

    :param model: a {@link #libsbml.Model Model} object.
    :param species_id2chebi_id: a mapping between species identifiers (string) and their ChEBI identifiers (string).
    :param ontology: ChEBI ontology.
    :param threshold: (Optional) A minimal number of reactions a species should participate in to become a ubiquitous one.
    The default value is {@link #UBIQUITOUS_THRESHOLD UBIQUITOUS_THRESHOLD}.
    :return: A set of ubiquitous species identifiers.
    """
    reactions = model.getListOfReactions()
    chebi2vote = {}
    for reaction in reactions:
        participants = get_metabolites(reaction)
        for element in participants:
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

    ubiquitous_chebi = set()
    for element, vote in chebi2vote.iteritems():
        if vote > threshold:
            ubiquitous_chebi.add(element[0])

    ubiquitous_chebi |= COMMON_UB_IDS
    return add_equivalent_chebi_ids(ontology, ubiquitous_chebi)


def get_ubiquitous_species_ids(model, ubiquitous_threshold=UBIQUITOUS_THRESHOLD):
    onto = parse(get_chebi())
    species_id2chebi_id = get_species_to_chebi(model, onto)
    cofactors = get_cofactors(onto)
    ubiquitous_chebi_ids = cofactors | get_ubiquitous_species_set(model, species_id2chebi_id, onto, ubiquitous_threshold)
    return {s.id for s in model.getListOfSpecies() if
            s.id in species_id2chebi_id and species_id2chebi_id[s.id] in ubiquitous_chebi_ids}


def get_cofactor_m_ids(model):
    chebi = parse(get_chebi())
    cofactors = get_cofactors(chebi) | COMMON_UB_IDS
    add_equivalent_chebi_ids(chebi, cofactors)
    s_id2ch_id = get_species_to_chebi(model, chebi, guess=True)
    return {s_id for s_id in s_id2ch_id.iterkeys() if s_id2ch_id[s_id] in cofactors}
