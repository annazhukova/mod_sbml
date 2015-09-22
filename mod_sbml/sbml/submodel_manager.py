from collections import Counter
from mod_sbml.sbml.sbml_manager import get_metabolites, get_reactants, get_products, get_modifiers, generate_unique_id
from mod_sbml.sbml.ubiquitous_manager import get_ubiquitous_chebi_ids, \
    select_metabolite_ids_by_term_ids

__author__ = 'anna'


def submodel(r_ids_to_keep, model):
    for r_id in [r.id for r in model.getListOfReactions() if r.id not in r_ids_to_keep]:
        model.removeReaction(r_id)
    remove_unused_species(model)
    remove_unused_compartments(model)


def remove_unused_species(model):
    s_ids_to_keep = set()
    for r in model.getListOfReactions():
        s_ids_to_keep |= get_metabolites(r, include_modifiers=True)
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
    biomass_r_ids = get_biomass_r_ids(model)
    r_ids_to_keep = {reaction.getId() for reaction in model.getListOfReactions()
                     if reaction.getId() not in biomass_r_ids}
    return submodel(r_ids_to_keep, model)


def get_biomass_r_ids(model):
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


def create_lumped_reaction(r_id2coeff, model, id_prefix='rl_', zero_threshold=1e-3, replace=True):
    new_r_id = generate_unique_id(model, id_=id_prefix)
    m_id2stoichiometry = compress_reaction_participants(model, r_id2coeff, zero_threshold=zero_threshold)
    new_r = model.createReaction()
    new_r.setId(new_r_id)
    new_r.setName('lumped reaction:%s' % ''.join(' %+g %s' % (coeff, r_id) for (r_id, coeff) in r_id2coeff.iteritems()))
    new_r.setReversible(False)
    for m_id, st in m_id2stoichiometry.iteritems():
        if st < 0:
            sr = new_r.createReactant()
            sr.setSpecies(m_id)
            sr.setStoichiometry(-st)
        elif st > 0:
            sr = new_r.createProduct()
            sr.setSpecies(m_id)
            sr.setStoichiometry(st)
    if replace:
        for r_id in r_id2coeff.iterkeys():
            model.removeReaction(r_id)
    return new_r_id


def compress_reaction_participants(model, r_id2coeff, zero_threshold=1e-3):
    m_id2stoichiometry = Counter()
    for r_id, coeff in r_id2coeff.iteritems():
        r = model.getReaction(r_id)
        m_id2stoichiometry.update({m_id: -st * coeff for (m_id, st) in get_reactants(r, stoichiometry=True)})
        m_id2stoichiometry.update({m_id: st * coeff for (m_id, st) in get_products(r, stoichiometry=True)})
    return {m_id: st for (m_id, st) in m_id2stoichiometry.iteritems() if abs(st) > zero_threshold}