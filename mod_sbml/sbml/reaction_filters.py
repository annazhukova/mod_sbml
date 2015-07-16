from mod_sbml.sbml.sbml_manager import get_metabolites, get_r_comps, get_reactants, get_products, get_modifiers, \
    get_genes, get_pathway_expression

__author__ = 'anna'


# by genes
def matches_genes(gene_collection, reaction):
    genes = get_genes(reaction)
    return set(genes) & set(gene_collection)


# by pathway
def matches_pathway(pathway_name, reaction):
    if not pathway_name:
        return False
    pathway_name = pathway_name.lower()
    pathways = get_pathway_expression(reaction)
    for pathway in pathways:
        if pathway and pathway.lower().find(pathway_name) != -1:
            return True
    return False


# by reaction attributes
def matches_ids(id_collection, reaction):
    return reaction.getId() in id_collection


def matches_name(name, reaction):
    if not name:
        return False
    name = name.lower()
    r_name = reaction.getName()
    return r_name and r_name.lower().find(name) != -1


# by compartment
def matches_compartment_id_weakly(compartment_id, reaction, model):
    return compartment_id in get_r_comps(reaction.getId(), model)


def matches_compartment_id(compartment_ids, reaction, model):
    return not (get_r_comps(reaction.getId(), model) - compartment_ids)


def matches_compartment_name_weakly(comp_name, reaction, model):
    if not comp_name:
        return False
    comp_name = comp_name.lower()
    for c_id in get_r_comps(reaction.getId(), model):
        c_name = model.getCompartment(c_id).getName()
        if c_name and c_name.lower().find(comp_name) != -1:
            return True
    return False


def matches_compartment_name(comp_name, reaction, model):
    if not comp_name:
        return False
    comp_name = comp_name.lower()
    c_ids = get_r_comps(reaction.getId(), model)
    for c_id in c_ids:
        c_name = model.getCompartment(c_id).getName()
        if not c_name or c_name.lower().find(comp_name) == -1:
            return False
    return len(c_ids) > 0


def is_not_transport(reaction, model):
    c_id = None
    participants = get_metabolites(reaction)
    for speciesId in participants:
        species = model.getSpecies(speciesId)
        compartment_id = species.getCompartment()
        if not compartment_id:
            return False
        if not c_id:
            c_id = compartment_id
        if compartment_id != c_id:
            return False
    return True


# by species
def matches_species_id(species_ids, reaction):
    return set(species_ids) & get_metabolites(reaction, include_modifiers=True)


def matches_species_name(name, reaction, model):
    if not name:
        return False
    name = name.lower()
    for speciesId in get_metabolites(reaction, include_modifiers=True):
        species = model.getSpecies(speciesId)
        if not species:
            continue
        species_name = species.getName()
        if species_name and species_name.lower().find(name) != -1:
            return True
    return False


def matches_reactant_id(s_id, reaction):
    return s_id in set(get_reactants(reaction))


def matches_product_id(s_id, reaction):
    return s_id in set(get_products(reaction))


def matches_modifier_id(s_id, reaction):
    return s_id in set(get_modifiers(reaction))


def matches_reactant_product_pair(reaction, s_id1, s_id2):
    reactants = set(get_reactants(reaction))
    products = set(get_products(reaction))
    if s_id1 in reactants:
        return s_id2 in products
    elif s_id2 in reactants:
        return s_id1 in products
    return False
