from collections import Counter

from chebi_annotator import get_species_to_chebi, CONJUGATE_ACID_BASE_RELATIONSHIPS, get_chebi, get_cofactors
from onto.obo_ontology import parse
from sbml.sbml_manager import get_metabolites


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
                 'chebi:73342', 'chebi:35780', 'chebi:26078'}

EQUIVALENT_TERM_RELATIONSHIPS = CONJUGATE_ACID_BASE_RELATIONSHIPS | {'is_tautomer_of'}


# # The function returns a set of identifiers of ubiquitous species participating in given reactions.
# The species in the model are divided into two groups: ubiquitous ones and the others.
# Ubiquitous species are those participating in more than {@link #threshold threshold number} of reactions.
# @param reactions A collection of {@link #libsbml.Reaction Reaction} objects.
# @param threshold (Optional) A minimal number of reactions a species should participate in to become a ubiquitous one.
# The default value is {@link #UBIQUITOUS_THRESHOLD UBIQUITOUS_THRESHOLD}.
# @return A set of ubiquitous species identifiers.
def get_ubiquitous_chebi_ids(model, species_id2chebi_id, ontology, threshold=UBIQUITOUS_THRESHOLD):
    chebi2vote = Counter()
    for reaction in model.getListOfReactions():
        for s_id in get_metabolites(reaction):
            # if we do not have a ChEBI annotation for it,
            # it will be considered ubiquitous anyway
            if s_id not in species_id2chebi_id:
                continue
            chebi_id = species_id2chebi_id[s_id]
            compartment = model.getSpecies(s_id).getCompartment()
            chebi2vote.update([(chebi_id, compartment)])

    ubiquitous_chebi = set()
    for c_id, vote in chebi2vote.iteritems():
        if vote > threshold:
            ubiquitous_chebi.add(c_id[0])
    ubiquitous_chebi |= COMMON_UB_IDS

    ubiquitous_chebi_new = set()
    for u_term_id in ubiquitous_chebi:
        u_term = ontology.get_term(u_term_id)
        if u_term:
            ubiquitous_chebi_new.add(u_term_id)
            ubiquitous_chebi_new |= {it.get_id() for it in
                                     ontology.get_equivalents(u_term, relationships=EQUIVALENT_TERM_RELATIONSHIPS)}
    return ubiquitous_chebi_new


def get_ubiquitous_species_ids(model, ubiquitous_threshold=UBIQUITOUS_THRESHOLD):
    onto = parse(get_chebi())
    species_id2chebi_id = get_species_to_chebi(model, onto)
    cofactors = get_cofactors(onto)
    ubiquitous_chebi_ids = cofactors | get_ubiquitous_chebi_ids(model, species_id2chebi_id, onto, ubiquitous_threshold)
    return {s.id for s in model.getListOfSpecies() if
            s.id in species_id2chebi_id and species_id2chebi_id[s.id] in ubiquitous_chebi_ids}