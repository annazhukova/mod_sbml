from collections import Counter

from mod_sbml.sbml.sbml_manager import get_metabolites, get_reactants, get_products, get_modifiers, \
    get_gene_association, set_gene_association

__author__ = 'anna'


def submodel_by_genes(model, genes_of_interest, keep_spontaneous_reactions=True):
    """
    Creates a submodel based on genes of interest. Reaction are kept
    if their gene associations can be satisfied only with the genes_of_interest.
    :param keep_spontaneous_reactions: whether the reactions that have no gene associations should be kept
    :param model: libsbml.Model
    :param genes_of_interest: a collection of genes of interest (must be in the same notation as in the model)
    :return:
    """
    r_ids = []
    spontaneous_r_ids = []
    for r in model.getListOfReactions():
        ga = get_gene_association(r, flatten=True, allowed_genes=genes_of_interest)
        if ga:
            r_ids.append(r.getId())
            set_gene_association(r, ga)
        elif ga is not None and keep_spontaneous_reactions:
            spontaneous_r_ids.append(r.getId())

    if spontaneous_r_ids:
        updated = True
        while updated:
            updated = False
            s_ids = set()
            for r in (model.getReaction(r_id) for r_id in r_ids):
                s_ids |= get_metabolites(r)

            for r in (model.getReaction(r_id) for r_id in spontaneous_r_ids):
                rs, ps = set(get_reactants(r)), set(get_products(r))
                if rs and rs == rs & s_ids or ps and ps == ps & s_ids:
                    r_ids.append(r.getId())
                    n = len(s_ids)
                    s_ids |= rs | ps
                    if len(s_ids) > n:
                        updated = True

    submodel(r_ids, model)


def submodel(r_ids_to_keep, model):
    for r_id in [r.id for r in model.getListOfReactions() if r.id not in r_ids_to_keep]:
        model.removeReaction(r_id)
    remove_unused_species(model)
    remove_unused_compartments(model)


def remove_unused_species(model, keep_modifiers=False):
    s_ids_to_keep = set()
    for r in model.getListOfReactions():
        s_ids_to_keep |= get_metabolites(r, include_modifiers=keep_modifiers)
    for s_id in [s.id for s in model.getListOfSpecies() if s.id not in s_ids_to_keep]:
        model.removeSpecies(s_id)


def remove_unused_compartments(model):
    c_ids_to_keep = {s.getCompartment() for s in model.getListOfSpecies()}
    c_ids_to_add = set()
    for c_id in c_ids_to_keep:
        while c_id:
            c_id = model.getCompartment(c_id).getOutside()
            if c_id and c_id not in c_ids_to_keep | c_ids_to_add:
                c_ids_to_add.add(c_id)
    c_ids_to_keep |= c_ids_to_add
    for c_id in [c.id for c in model.getListOfCompartments() if c.id not in c_ids_to_keep]:
        model.removeCompartment(c_id)


def biomassless_model(model):
    biomass_r_ids = get_biomass_r_ids(model)
    r_ids_to_keep = {reaction.getId() for reaction in model.getListOfReactions()
                     if reaction.getId() not in biomass_r_ids}
    return submodel(r_ids_to_keep, model)


def get_biomass_r_ids(model):
    """
    Find biomass-related reactions in the model.
    :param model: libsbml.Model model of interest
    :return: set of biomass-related reaction ids
    """
    has_biomass = lambda s: s and (-1 != s.lower().find('biomass') or -1 != s.lower().find('growth'))
    biomass_s_ids = {species.getId() for species in model.getListOfSpecies()
                     if has_biomass(species.getName()) or has_biomass(species.getId())}
    return {reaction.getId() for reaction in model.getListOfReactions()
            if has_biomass(reaction.getName()) or has_biomass(reaction.getId())
            or biomass_s_ids & get_metabolites(reaction, include_modifiers=False)}


def remove_species(model, s_ids_to_remove):
    for r in model.getListOfReactions():
        for iterator, remove in ((get_reactants(r), r.removeReactant), (get_products(r), r.removeProduct),
                                 (get_modifiers(r), r.removeModifier)):
            for s_id in s_ids_to_remove & set(iterator):
                remove(s_id)
        if not r.getNumProducts() and not r.getNumReactants():
            model.removeReaction(r.getId())
    for s_id in s_ids_to_remove:
        model.removeSpecies(s_id)
    remove_unused_compartments(model)


def compress_reaction_participants(model, r_id2coeff, zero_threshold=1e-3):
    m_id2stoichiometry = Counter()
    for r_id, coeff in r_id2coeff.items():
        r = model.getReaction(r_id)
        m_id2stoichiometry.update({m_id: -st * coeff for (m_id, st) in get_reactants(r, stoichiometry=True)})
        m_id2stoichiometry.update({m_id: st * coeff for (m_id, st) in get_products(r, stoichiometry=True)})
    return {m_id: -st for (m_id, st) in m_id2stoichiometry.items() if st < -zero_threshold},\
           {m_id: st for (m_id, st) in m_id2stoichiometry.items() if st > zero_threshold}
