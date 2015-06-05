from collections import defaultdict
import logging

import libsbml

from onto.obo_ontology import to_identifiers_org_format


SBO_MATERIAL_ENTITY = "SBO:0000240"

__author__ = 'anna'

URN_MIRIAM = "urn:miriam:"

IDENTIFIERS_ORG = "http://identifiers.org/"

EC_PREFIX = "EC Number:"

GA_PREFIX = "GENE_ASSOCIATION:"

PATHWAY_PREFIX = "SUBSYSTEM:"

FORMULA_PREFIX = "FORMULA:"


def get_annotation_term_of_type(element, qualifier_type):
    cv_terms = element.getCVTerms()
    if cv_terms:
        for i in xrange(cv_terms.getSize()):
            term = cv_terms.get(i)
            if libsbml.BIOLOGICAL_QUALIFIER == term.getQualifierType() and qualifier_type == term.getBiologicalQualifierType():
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
        term = libsbml.CVTerm(libsbml.BIOLOGICAL_QUALIFIER)
        term.setBiologicalQualifierType(qualifier)
        term.addResource(annotation)
        element.addCVTerm(term, True)
    else:
        term.addResource(annotation)


def remove_annotation(element, qualifier, annotation):
    if not element.isSetMetaId():
        element.setMetaId("m_{0}".format(element.getId()))
    for term in get_annotation_term_of_type(element, qualifier):
        term.removeResource(annotation)


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


def get_subsystem2r_ids(sbml=None, model=None):
    subsystem2r_ids = defaultdict(set)
    no_pathway_r_ids = set()
    if not model and not sbml:
        raise ValueError("Either sbml or model parameter should be specified")
    if not model:
        input_doc = libsbml.SBMLReader().readSBML(sbml)
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


def get_formulas(species):
    result = set()
    node = species.getNotes()
    _get_prefixed_notes_value(node, result, FORMULA_PREFIX)
    return result


def get_subsystem(reaction):
    result = set()
    node = reaction.getNotes()
    _get_prefixed_notes_value(node, result, PATHWAY_PREFIX)
    return result


def get_pathway2r_ids(sbml=None, model=None):
    pw2r_ids = defaultdict(set)
    no_pw_r_ids = set()
    if not model and not sbml:
        raise ValueError("Either sbml or model parameter should be specified")
    if not model:
        input_doc = libsbml.SBMLReader().readSBML(sbml)
        model = input_doc.getModel()
    for r in model.getListOfReactions():
        found = False
        for annotation in get_annotations(r, libsbml.BQB_IS_PART_OF):
            if annotation.find("kegg.pathway") != -1:
                pw2r_ids[annotation.replace("kegg.pathway:", '')].add(r.getId())
                found = True
        if not found:
            no_pw_r_ids.add(r.getId())
    # Should find the following pathways:
    # 1. citrate cycle (TCA cycle) (hsa00020);
    # 2. fatty acid metabolism (hsa01212) and 3 sub-pathways:
    #   (a) fatty acid biosynthesis (hsa00061);
    #   (b) fatty acid elongation (hsa00062);
    #   (c) fatty acid degradation (hsa00071);
    # 3. valine, leucine and isoleucine degradation (hsa00280);
    # 4. synthesis and degradation of ketone bodies (hsa00072).

    # Let's remove the sub-pathways
    if "path:hsa01212" in pw2r_ids:
        for key in ("path:hsa00061", "path:hsa00062", "path:hsa00071"):
            if key in pw2r_ids:
                del pw2r_ids[key]

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
    for annotation in get_annotations(r, libsbml.BQB_IS):
        if annotation.find("kegg.reaction") != -1:
            return annotation.replace("kegg.reaction:", '')
    return None


def get_kegg_m_id(m):
    for annotation in get_annotations(m, libsbml.BQB_IS):
        if annotation.find("kegg.compound") != -1:
            return annotation.replace("kegg.compound:", '')
    return None


def get_chebi_id(m):
    for annotation in get_annotations(m, libsbml.BQB_IS):
        if annotation.lower().find("chebi") != -1:
            return annotation.lower().replace("obo.chebi:", '').replace("chebi:chebi:", 'chebi:')
    return None


def get_taxonomy(model):
    occurs_in = get_qualifier_values(model.getAnnotation(), libsbml.BQB_OCCURS_IN)
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


def get_metabolites(reaction, stoichiometry=False):
    return set(get_reactants(reaction, stoichiometry)) | set(get_products(reaction, stoichiometry))


def get_r_comps(r_id, model):
    r = model.getReaction(r_id)
    return {model.getSpecies(s_id).getCompartment() for s_id in get_metabolites(r)}


def create_species(model, compartment_id, name=None, bound=False, id_=None):
    new_species = model.createSpecies()
    id_ = generate_unique_id(model, id_ if id_ else "s")
    if libsbml.LIBSBML_OPERATION_SUCCESS != new_species.setId(id_):
        logging.error("species  %s creation error" % id_)
    if name:
        new_species.setName(name)
    new_species.setCompartment(compartment_id)
    new_species.setSBOTerm(SBO_MATERIAL_ENTITY)
    new_species.setBoundaryCondition(bound)
    return new_species


def create_compartment(model, name, outside=None, term_id=None, id_=None):
    new_comp = model.createCompartment()
    id_ = generate_unique_id(model, id_ if id_ else "c")
    new_comp.setId(id_)
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


def generate_unique_id(model, id_=None):
    if not id_:
        id_ = "new_"
    else:
        id_ = ''.join(e for e in id_ if e.isalnum())
        if id_[0].isdigit():
            id_ = 's_' + id_
        if not model.getElementBySId(id_):
            return id_
    i = 0
    while model.getElementBySId("%s_%d" % (id_, i)):
        i += 1
    return "%s_%d" % (id_, i)


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


def create_reaction(model, rs, ps, name=None, reversible=True, _id=None):
    def add_r_elements(elements, creator):
        for (m_id, st) in elements:
            element = creator()
            element.setSpecies(m_id)
            element.setStoichiometry(st)

    r = model.createReaction()
    r_id = generate_unique_id(model, _id if _id else 'r')
    r.setName(name)
    r.setReversible(reversible)
    if libsbml.LIBSBML_OPERATION_SUCCESS != r.setId(r_id):
        logging.error("reaction  ", r_id, " creation error")
    add_r_elements(rs, r.createReactant)
    add_r_elements(ps, r.createProduct)
    return r


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