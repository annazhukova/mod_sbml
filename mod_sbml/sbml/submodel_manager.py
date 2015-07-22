from mod_sbml.sbml.sbml_manager import get_metabolites, get_reactants, get_products
from mod_sbml.sbml.ubiquitous_manager import get_ubiquitous_chebi_ids, \
    select_metabolite_ids_by_term_ids

__author__ = 'anna'


def submodel(r_ids_to_keep, model):
    for r_id in [r.id for r in model.getListOfReactions() if r.id not in r_ids_to_keep]:
        model.removeReaction(r_id)
    s_ids_to_keep = set()
    for r in model.getListOfReactions():
        s_ids_to_keep |= get_metabolites(r, include_modifiers=True)
    for s_id in [s.id for s in model.getListOfSpecies() if s.id not in s_ids_to_keep]:
        model.removeSpecies(s_id)
    remove_unused_compartments(model)


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


def specific_metabolite_model(model, cofactor_m_ids=None):
    if not cofactor_m_ids:
        ub_ch_ids = get_ubiquitous_chebi_ids(add_common=True, add_cofactors=True)
        cofactor_m_ids = select_metabolite_ids_by_term_ids(model, ub_ch_ids)
    r_ids_to_remove = []
    for r in model.getListOfReactions():
        for m_id in set(get_reactants(r)) & cofactor_m_ids:
            r.removeReactant(m_id)
        for m_id in set(get_products(r)) & cofactor_m_ids:
            r.removeProduct(m_id)
        if not r.getNumProducts() and not r.getNumReactants():
            r_ids_to_remove.append(r.id)
    for r_id in r_ids_to_remove:
        model.removeReaction(r_id)
    for m_id in cofactor_m_ids:
        model.removeSpecies(m_id)
    remove_unused_compartments(model)


def biomassless_model(model):
    has_biomass = lambda s: s and -1 != s.lower().find('biomass')
    biomass_s_ids = {species.getId() for species in model.getListOfSpecies()
                     if has_biomass(species.getName()) or has_biomass(species.getId())}
    r_ids_to_keep = {reaction.getId() for reaction in model.getListOfReactions()
                     if not has_biomass(reaction.getName()) and not has_biomass(reaction.getId())
                     and not biomass_s_ids & get_metabolites(reaction, include_modifiers=True)}
    return submodel(r_ids_to_keep, model)

