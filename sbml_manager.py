from collections import defaultdict
import logging
from libsbml import *
from pathway_manager import get_relevant_pathway_info

SBO_MATERIAL_ENTITY = "SBO:0000240"

__author__ = 'anna'

URN_MIRIAM = "urn:miriam:"

IDENTIFIERS_ORG = "http://identifiers.org/"

EC_PREFIX = "EC Number:"

GA_PREFIX = "GENE_ASSOCIATION:"

PATHWAY_PREFIX = "SUBSYSTEM:"


def get_annotation_term_of_type(element, qualifier_type):
    cv_terms = element.getCVTerms()
    if cv_terms:
        for i in xrange(cv_terms.getSize()):
            term = cv_terms.get(i)
            if BIOLOGICAL_QUALIFIER == term.getQualifierType() and qualifier_type == term.getBiologicalQualifierType():
                yield term


def get_qualifier_values(element, qualifier_type):
    for term in get_annotation_term_of_type(element, qualifier_type):
        for i in xrange(term.getNumResources()):
            yield term.getResourceURI(i).replace("%3A", ":")


def add_annotation(element, qualifier, annotation):
    if not element.isSetMetaId():
        element.setMetaId("m_{0}".format(element.getId()))
    term = next(get_annotation_term_of_type(element, qualifier), None)
    if not term:
        term = CVTerm(BIOLOGICAL_QUALIFIER)
        term.setBiologicalQualifierType(qualifier)
        term.addResource(annotation)
        element.addCVTerm(term, True)
    else:
        term.addResource(annotation)


def miriam_to_term_id(urn):
    urn = urn.strip()
    miriam_prefix = URN_MIRIAM
    identifiers_org_prefix = IDENTIFIERS_ORG
    if 0 == urn.find(miriam_prefix):
        urn = urn[len(miriam_prefix):]
    elif 0 == urn.find(identifiers_org_prefix):
        urn = urn[len(identifiers_org_prefix):]
    return urn.strip().replace("/", ":")


def get_annotations(entity, qualifier):
    return (miriam_to_term_id(it) for it in get_qualifier_values(entity, qualifier))


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


def get_subsystem2r_ids(sbml):
    subsystem2r_ids = defaultdict(set)
    no_pathway_r_ids = set()
    input_doc = SBMLReader().readSBML(sbml)
    model = input_doc.getModel()
    for r in model.getListOfReactions():
        result = set()
        node = r.getNotes()
        _get_prefixed_notes_value(node, result, PATHWAY_PREFIX)
        for pw in result:
            subsystem2r_ids[pw].add(r.getId())
        if not result:
            no_pathway_r_ids.add(r.getId())
    return subsystem2r_ids, no_pathway_r_ids


def get_pathway2r_ids(sbml):
    pw2r_ids = defaultdict(set)
    no_pw_r_ids = set()
    input_doc = SBMLReader().readSBML(sbml)
    model = input_doc.getModel()
    for r in model.getListOfReactions():
        found = False
        for annotation in get_annotations(r, BQB_IS_PART_OF):
            if annotation.find("kegg.pathway") != -1:
                pw2r_ids[annotation.replace("kegg.pathway:", '')].add(r.getId())
                found = True
        if not found:
            no_pw_r_ids.add(r.getId())
    return pw2r_ids, no_pw_r_ids


def get_kegg_r_id2r_ids(model):
    kegg_r_id2r_ids = defaultdict(set)
    for r in model.getListOfReactions():
        k_r_id = get_kegg_r_id(r)
        if k_r_id:
            kegg_r_id2r_ids[k_r_id].add(r.getId())
    return kegg_r_id2r_ids


def get_kegg_m_id2m_ids(model):
    kegg_m_id2m_ids = defaultdict(set)
    for m in model.getListOfSpecies():
        k_m_id = get_kegg_m_id(m)
        if k_m_id:
            kegg_m_id2m_ids[k_m_id].add(m.getId())
    return kegg_m_id2m_ids


def get_kegg_r_id(r):
    for annotation in get_annotations(r, BQB_IS):
        if annotation.find("kegg.reaction") != -1:
            return annotation.replace("kegg.reaction:", '')
    return None


def get_kegg_m_id(m):
    for annotation in get_annotations(m, BQB_IS):
        if annotation.find("kegg.compound") != -1:
            return annotation.replace("kegg.compound:", '')
    return None


def get_chebi_id(m):
    for annotation in get_annotations(m, BQB_IS):
        if annotation.lower().find("chebi") != -1:
            return annotation.lower().replace("obo.chebi:", '').replace("chebi:chebi:", 'chebi:')
    return None


def copy_sbml(sbml_in, sbml_out):
    input_doc = SBMLReader().readSBML(sbml_in)
    SBMLWriter().writeSBMLToFile(input_doc, sbml_out)


def annotate_with_pathways(org, sbml, kegg_r_ids, kegg_r_id2r_ids, threshold=0.5):
    pw2name_rs_ratio = get_relevant_pathway_info(org, {"rn:" + r_id for r_id in kegg_r_ids}, threshold=threshold)
    input_doc = SBMLReader().readSBML(sbml)
    model = input_doc.getModel()
    for pw, (name, rns, ratio) in pw2name_rs_ratio.iteritems():
        pw = pw.replace("pathway:", "")
        for kegg_r_id in rns:
            kegg_r_id = kegg_r_id.replace("rn:", "")
            for r_id in kegg_r_id2r_ids(kegg_r_id):
                r = model.getElementBySId(r_id)
                if r:
                    add_annotation(r, BQB_IS_PART_OF, "http://identifiers.org/kegg.pathway/%s" % pw)
                    r.appendNotes("<html:body><html:p>SUBSYSTEM: {1} ({0})</html:p></html:body>".format(pw, name))
    SBMLWriter().writeSBMLToFile(input_doc, sbml)


def get_stoichiometry(species_ref):
    result = species_ref.getStoichiometry()
    if not result:
        result = species_ref.getStoichiometryMath()
    if not result:
        return 1
    return result


def get_reactants(reaction, stoichiometry=False):
    if stoichiometry:
        return {(species_ref.getSpecies(), get_stoichiometry(species_ref)) for species_ref in
                reaction.getListOfReactants()}
    else:
        return {species_ref.getSpecies() for species_ref in reaction.getListOfReactants()}


def get_products(reaction, stoichiometry=False):
    if stoichiometry:
        return {(species_ref.getSpecies(), get_stoichiometry(species_ref)) for species_ref in
                reaction.getListOfProducts()}
    else:
        return {species_ref.getSpecies() for species_ref in reaction.getListOfProducts()}


def get_metabolites(reaction, stoichiometry=False):
    return get_reactants(reaction, stoichiometry) | get_products(reaction, stoichiometry)


def get_r_comp(r_id, model):
    r = model.getReaction(r_id)
    return {model.getSpecies(s_id).getCompartment() for s_id in get_metabolites(r)}


def create_species(model, compartment_id, name=None, bound=False, id_=None):
    new_species = model.createSpecies()
    id_ = generate_unique_id(model, id_)
    if LIBSBML_OPERATION_SUCCESS != new_species.setId(id_):
        logging.error("species  %s creation error" % id_)
    new_species.setName(name)
    new_species.setCompartment(compartment_id)
    new_species.setSBOTerm(SBO_MATERIAL_ENTITY)
    new_species.setBoundaryCondition(bound)
    return new_species


def generate_unique_id(model, id_=None):
    if not id_:
        id_ = "new_"
    else:
        id_ = ''.join(e for e in id_ if e.isalnum())
        if not model.getElementBySId(id_):
            return id_
    i = 0
    while model.getElementBySId("%s_%d" % (id_, i)):
        i += 1
    return "%s_%d" % (id_, i)


def get_reaction_ids(model, s_id):
    return {r.id for r in model.getListOfReactions() if s_id in get_metabolites(r)}


def get_reactions_by_comp(model, comps, strict=False):
    r_ids = set()
    for r in model.getListOfReactions():
        s_ids = get_metabolites(r)
        if strict and not {model.getSpecies(s_id).getCompartment() for s_id in s_ids} - comps or \
                        not strict and {model.getSpecies(s_id).getCompartment() for s_id in s_ids} & comps:
            r_ids.add(r.id)
    return r_ids


def submodel(r_ids_to_keep, model):
    for r_id in [r.id for r in model.getListOfReactions() if not r.id in r_ids_to_keep]:
        model.removeReaction(r_id)
    s_ids_to_keep = set()
    for r in model.getListOfReactions():
        s_ids_to_keep |= get_metabolites(r)
    for s_id in [s.id for s in model.getListOfSpecies() if not s.id in s_ids_to_keep]:
        model.removeSpecies(s_id)
    c_ids_to_keep = {s.getCompartment() for s in model.getListOfSpecies()}
    c_ids_to_add = set()
    for c_id in c_ids_to_keep:
        while c_id:
            c_id = model.getCompartment(c_id).getOutside()
            if c_id and not c_id in c_ids_to_keep | c_ids_to_add:
                c_ids_to_add.add(c_id)
    c_ids_to_keep |= c_ids_to_add
    for c_id in [c.id for c in model.getListOfCompartments() if not c.id in c_ids_to_keep]:
        model.removeCompartment(c_id)