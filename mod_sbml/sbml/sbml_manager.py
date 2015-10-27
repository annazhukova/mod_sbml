from collections import defaultdict
import logging
import os

import libsbml

from mod_sbml.annotation.miriam_converter import to_identifiers_org_format
from mod_sbml.annotation.rdf_annotation_helper import get_qualifier_values, add_annotation

SBO_MATERIAL_ENTITY = "SBO:0000240"

__author__ = 'anna'

EC_PREFIX = "EC Number:"

GA_PREFIX = "GENE_ASSOCIATION:"

PATHWAY_PREFIX = "SUBSYSTEM:"

FORMULA_PREFIX = "FORMULA:"


def get_model_name(sbml):
    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()
    model_name = model.getName() if model.getName() else model.getId()
    if not model_name or not model_name.strip():
        return os.path.splitext(os.path.basename(sbml))[0]
    return model_name.strip()


def _get_prefixed_notes_value(notes, result, prefix):
    if not notes:
        return
    for i in xrange(0, notes.getNumChildren()):
        child = notes.getChild(i)
        note = child.getCharacters()
        if note:
            start = note.find(prefix)
            if start != -1:
                start += len(prefix)
                result.add(note[start:len(note)].strip())
        _get_prefixed_notes_value(child, result, prefix)


def get_ec_numbers(reaction):
    result = set()
    node = reaction.getNotes()
    _get_prefixed_notes_value(node, result, EC_PREFIX)
    return result


def get_gene_association(reaction):
    result = set()
    node = reaction.getNotes()
    _get_prefixed_notes_value(node, result, GA_PREFIX)
    if result:
        return result.pop()
    return ''


def get_genes(reaction):
    return gene_association2genes(get_gene_association(reaction))


def gene_association2genes(gene_association):
    genes = []
    if gene_association:
        for g0 in gene_association.split('('):
            for g1 in g0.split(')'):
                for g2 in g1.split('and'):
                    for g3 in g2.split('or'):
                        for g4 in g3.split('xor'):
                            for g5 in g4.split('not'):
                                g5 = g5.replace(' ', '')
                                if g5:
                                    genes.append(g5)
    genes.sort()
    return genes


def get_pathway_expression(reaction):
    result = set()
    node = reaction.getNotes()
    _get_prefixed_notes_value(node, result, PATHWAY_PREFIX)
    return result


def get_subsystem2r_ids(sbml=None, model=None):
    subsystem2r_ids = defaultdict(set)
    no_pathway_r_ids = set()
    if not model and not sbml:
        raise ValueError("Either sbml or model parameter should be specified")
    if not model:
        input_doc = libsbml.SBMLReader().readSBML(sbml)
        model = input_doc.getModel()
    for r in model.getListOfReactions():
        result = get_pathway_expression(r)
        for pw in result:
            subsystem2r_ids[pw].add(r.getId())
        if not result:
            no_pathway_r_ids.add(r.getId())
    return subsystem2r_ids, no_pathway_r_ids


def get_formulas(species):
    result = set()
    node = species.getNotes()
    _get_prefixed_notes_value(node, result, FORMULA_PREFIX)
    return {formula.strip() for formula in result if formula and '.' != formula.strip()}


def get_subsystem(reaction):
    result = set()
    node = reaction.getNotes()
    _get_prefixed_notes_value(node, result, PATHWAY_PREFIX)
    return result


def get_taxonomy(model):
    occurs_in = get_qualifier_values(model, libsbml.BQB_OCCURS_IN)
    for it in occurs_in:
        start = it.find("taxonomy")
        if start != -1:
            return it[start + len("taxonomy:"):].strip()
    return None


def copy_sbml(sbml_in, sbml_out):
    input_doc = libsbml.SBMLReader().readSBML(sbml_in)
    libsbml.SBMLWriter().writeSBMLToFile(input_doc, sbml_out)


def get_stoichiometry(species_ref):
    result = species_ref.getStoichiometry()
    if result:
        return result
    result = species_ref.getStoichiometryMath()
    if result:
        return result
    return 1


def get_reactants(reaction, stoichiometry=False):
    if stoichiometry:
        return ((species_ref.getSpecies(), get_stoichiometry(species_ref)) for species_ref in
                reaction.getListOfReactants())
    else:
        return (species_ref.getSpecies() for species_ref in reaction.getListOfReactants())


def get_products(reaction, stoichiometry=False):
    if stoichiometry:
        return ((species_ref.getSpecies(), get_stoichiometry(species_ref)) for species_ref in
                reaction.getListOfProducts())
    else:
        return (species_ref.getSpecies() for species_ref in reaction.getListOfProducts())


def get_modifiers(reaction, stoichiometry=False):
    if stoichiometry:
        return ((species_ref.getSpecies(), get_stoichiometry(species_ref)) for species_ref in
                reaction.getListOfModifiers())
    else:
        return (species_ref.getSpecies() for species_ref in reaction.getListOfModifiers())


def get_metabolites(reaction, stoichiometry=False, include_modifiers=False):
    result = set(get_reactants(reaction, stoichiometry)) | set(get_products(reaction, stoichiometry))
    if include_modifiers:
        result |= set(get_modifiers(reaction, stoichiometry))
    return result


def get_r_comps(r_id, model):
    r = model.getReaction(r_id)
    return {model.getSpecies(s_id).getCompartment() for s_id in get_metabolites(r, include_modifiers=True)}


def create_species(model, compartment_id, name=None, bound=False, id_=None, type_id=None, sbo_id=None, kegg_id=None,
                   chebi_id=None):
    new_species = model.createSpecies()
    id_ = generate_unique_id(model, id_ if id_ else "s")
    if libsbml.LIBSBML_OPERATION_SUCCESS != new_species.setId(id_):
        logging.error("species  %s creation error" % id_)
    if name:
        new_species.setName(name)
    new_species.setCompartment(compartment_id)
    if type_id:
        new_species.setSpeciesType(type_id)
    new_species.setSBOTerm(sbo_id if sbo_id else SBO_MATERIAL_ENTITY)
    new_species.setBoundaryCondition(bound)
    if kegg_id:
        add_annotation(new_species, libsbml.BQB_IS, "http://identifiers.org/kegg.compound/%s" % kegg_id.upper())
    if chebi_id:
        add_annotation(new_species, libsbml.BQB_IS, "http://identifiers.org/obo.chebi/%s" % chebi_id.upper())
    return new_species


def create_compartment(model, name=None, outside=None, term_id=None, id_=None):
    new_comp = model.createCompartment()
    id_ = generate_unique_id(model, id_ if id_ else "c")
    if libsbml.LIBSBML_OPERATION_SUCCESS != new_comp.setId(id_):
        logging.error("compartment %s creation error" % id_)
    if name:
        new_comp.setName(name)
    if outside:
        new_comp.setOutside(outside)
    if term_id:
        add_annotation(new_comp, libsbml.BQB_IS, to_identifiers_org_format(term_id, "obo.go"))
    return new_comp


def reverse_reaction(r):
    rs = list(r.getListOfReactants())
    ps = list(r.getListOfProducts())
    for m in rs:
        new_m = r.createProduct()
        new_m.setSpecies(m.getSpecies())
        new_m.setStoichiometry(m.getStoichiometry())
        r.removeReactant(m.getSpecies())
    for m in ps:
        new_m = r.createReactant()
        new_m.setSpecies(m.getSpecies())
        new_m.setStoichiometry(m.getStoichiometry())
        r.removeProduct(m.getSpecies())


def generate_unique_id(model, id_=None, i=0):
    if not id_:
        id_ = 's_'
    else:
        id_ = ''.join(e for e in id_ if e.isalnum() or '_' == e)
        if id_[0].isdigit() or '_' == id_[0]:
            id_ = 's_' + id_
    if not model.getElementBySId(id_):
        return id_
    while model.getElementBySId("%s%d" % (id_, i)):
        i += 1
    return "%s%d" % (id_, i)


def get_r_ids_by_s_ids(model, s_ids):
    return {r.id for r in model.getListOfReactions() if set(s_ids) & get_metabolites(r)}


def get_r_ids_by_comp(model, comps, strict=False):
    r_ids = set()
    for r in model.getListOfReactions():
        s_ids = get_metabolites(r)
        if strict and not {model.getSpecies(s_id).getCompartment() for s_id in s_ids} - comps or not strict and {
            model.getSpecies(s_id).getCompartment() for s_id in s_ids} & comps:
            r_ids.add(r.id)
    return r_ids


def get_pathway_by_species(s_ids, model, ubiquitous_s_ids, blocked_r_ids=None):
    s_id2r_ids, r_id2s_ids = defaultdict(set), {}
    for r in model.getListOfReactions():
        if blocked_r_ids and r.id in blocked_r_ids:
            continue
        sp_ids = get_metabolites(r) - ubiquitous_s_ids
        r_id2s_ids[r.id] = sp_ids
        for sp_id in sp_ids:
            s_id2r_ids[sp_id].add(r.id)
    s_ids_to_process = set(s_ids)
    s_ids_processed = set()
    r_ids = set()
    while s_ids_to_process:
        s_id = s_ids_to_process.pop()
        s_ids_processed.add(s_id)
        if s_id not in s_id2r_ids:
            continue
        r_ids_to_add = s_id2r_ids[s_id] - r_ids
        r_ids |= r_ids_to_add
        for r_id in r_ids_to_add:
            s_ids_to_process |= (r_id2s_ids[r_id] - s_ids_processed)

    for s in sorted((model.getSpecies(s_id) for s_id in s_ids_processed), key=lambda s: s.name):
        print(s.name, s.id, len(s_id2r_ids[s.id]))

    print("________________________________________")
    for s in sorted((model.getSpecies(s_id) for s_id in ubiquitous_s_ids), key=lambda s: s.name):
        print(s.name, s.id)

    return r_ids


def create_reaction(model, r_id2st, p_id2st, name=None, reversible=True, id_=None, kegg_id=None):
    new_r_id = generate_unique_id(model, id_=id_)
    new_r = model.createReaction()
    if libsbml.LIBSBML_OPERATION_SUCCESS != new_r.setId(new_r_id):
        logging.error("reaction %s creation error" % new_r_id)
    if name:
        new_r.setName(name)
    new_r.setReversible(reversible)
    for m_id, st in r_id2st.iteritems():
        sr = new_r.createReactant()
        sr.setSpecies(m_id)
        sr.setStoichiometry(st)
    for m_id, st in p_id2st.iteritems():
        sr = new_r.createProduct()
        sr.setSpecies(m_id)
        sr.setStoichiometry(st)
    if kegg_id:
        add_annotation(new_r, libsbml.BQB_IS, "http://identifiers.org/kegg.reaction/%s" % kegg_id.upper())
    return new_r


def find_reaction_by_reactant_product_names(model, r_name, p_name):
    r_name, p_name = r_name.lower(), p_name.lower()
    for r in model.getListOfReactions():
        r_found, p_found, rev = False, False, False
        for m in (model.getSpecies(s_id) for s_id in get_reactants(r)):
            if m.name.lower().find(r_name) != -1:
                r_found = True
                break
            if m.name.lower().find(p_name) != -1:
                r_found = True
                rev = True
                break
        if not r_found:
            continue
        for m in (model.getSpecies(s_id) for s_id in get_products(r)):
            if rev and m.name.lower().find(r_name) != -1:
                p_found = True
                break
            if not rev and m.name.lower().find(p_name) != -1:
                p_found = True
                break
        if p_found:
            return r
    return None
