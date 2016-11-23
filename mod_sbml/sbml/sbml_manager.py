from collections import defaultdict
import logging
import os

import libsbml
import pyparsing as pp

OR = ["OR", "or", "Or"]

AND = ["AND", "and", "And"]

SBO_MATERIAL_ENTITY = "SBO:0000240"

__author__ = 'anna'

EC_PREFIX = "EC Number:"

GA_PREFIX = "GENE_ASSOCIATION:"

PATHWAY_PREFIX = "SUBSYSTEM:"

FORMULA_PREFIX = "FORMULA:"


def get_model_name(sbml=None, model=None):
    if not model:
        doc = libsbml.SBMLReader().readSBML(sbml)
        model = doc.getModel()
    model_name = model.getName() if model.getName() else model.getId()
    if not model_name or not model_name.strip():
        return os.path.splitext(os.path.basename(sbml))[0]
    return model_name.replace('_', ' ').strip()


def _get_prefixed_notes_value(notes, result, prefix):
    if not notes:
        return
    for i in range(0, notes.getNumChildren()):
        child = notes.getChild(i)
        note = child.getCharacters()
        if note:
            start = note.find(prefix)
            if start != -1:
                start += len(prefix)
                result.add(note[start:].strip())
        _get_prefixed_notes_value(child, result, prefix)


def get_ec_numbers(reaction):
    result = set()
    node = reaction.getNotes()
    _get_prefixed_notes_value(node, result, EC_PREFIX)
    return result


def _remove_duplicates(expression, flatten=False):
    if not isinstance(expression, list):
        return expression
    if len(expression) % 2 == 0:
        raise ValueError("The expression was supposed to be in a form [i_1, op, i_2, op, ... op, i_n], "
                         "got even number of elements instead")
    genes = {_remove_duplicates(it, flatten=flatten) for it in expression[::2]}
    if len(genes) == 1:
        return genes.pop()

    operand = expression[1]
    result = []
    for gene in genes:
        result.append(gene)
        result.append(operand)
    return '(%s)' % ' '.join(result[: -1]) if flatten else tuple(result[: -1])


def _filter(expression, allowed_values, flatten=False):
    if not allowed_values:
        return expression
    if not isinstance(expression, tuple):
        return expression if expression in allowed_values else None
    if len(expression) % 2 == 0:
        raise ValueError("The expression was supposed to be in a form [i_1, op, i_2, op, ... op, i_n], "
                         "got even number of elements instead")
    genes = {_filter(it, allowed_values, flatten) for it in expression[::2]}

    operand = expression[1]
    filtered_genes = [gene for gene in genes if gene is not None]

    if operand in AND and len(filtered_genes) < len(genes) or not filtered_genes:
        return None

    if len(filtered_genes) == 1:
        return genes.pop()

    result = []
    for gene in filtered_genes:
        result.append(gene)
        result.append(operand)
    return '(%s)' % ' '.join(result[: -1]) if flatten else tuple(result[: -1])


def parse_gene_association(ga, gene_parse_action=None, flatten=True):
    if not ga:
        return ga

    gene = pp.Word(initChars=pp.alphanums + "_.-")

    if gene_parse_action:
        gene = gene.setParseAction(gene_parse_action)

    and_op = pp.oneOf(AND)
    or_op = pp.oneOf(OR)

    expr = pp.operatorPrecedence(gene, [
        (and_op, 2, pp.opAssoc.LEFT,),
        (or_op, 2, pp.opAssoc.LEFT,)
    ])

    res = expr.parseString(ga, parseAll=True).asList()
    while len(res) == 1:
        res = res[0]
    return _remove_duplicates(res, flatten=flatten)


def get_gene_association(reaction, gene_parse_action=None, flatten=True, allowed_genes=None):
    """
    Extracts gene association encoded as a note prefixed with GENE_ASSOCIATION:.
    :param allowed_genes: if allowed_genes are specified the resulting gene association will be filtered
    to only contain them
    :param reaction: the reaction of interest (libsbml.Reaction)
    :param gene_parse_action: if you need to do something with the genes, e.g. convert them to a different format,
    add it as a gene parse action. The gene parse action takes a term as an input (the gene value is term[0])
    and returns the new value as an output. For example, if you have a function entrez2symbol that
    converts ENTREZ gene ID to SYMBOL, specify gene_parse_action as lambda term: entrez2symbol(term[0]).
    If gene_parse_action is None (default), no modification to the genes is made.
    :param flatten: whether to return the gene association as a string (True), e.g. '((3906 and 2683) or 8704)',
    or as a list (False), e.g. [['3906', 'and', '2683'], 'or', '8704'].
    :return: the gene association as a string (if flatten is True), e.g. '((3906 and 2683) or 8704)',
    or as a list (if flatten is False), e.g. [['3906', 'and', '2683'], 'or', '8704']; empty for spontaneous reactions;
    or None if there was a parsing problem or if the needed genes are not present among the allowed ones.
    """
    result = set()
    node = reaction.getNotes()
    _get_prefixed_notes_value(node, result, GA_PREFIX)
    if result:
        ga = result.pop()
        try:
            ga = parse_gene_association(ga, gene_parse_action=gene_parse_action,
                                        flatten=False if allowed_genes else flatten)
            if ga and allowed_genes:
                f_ga = _filter(ga, allowed_genes, flatten=flatten)
                return f_ga if f_ga else None
            return ga
        except pp.ParseBaseException as e:
            logging.error('Ignoring the gene association for %s as it is malformed: %s' % (reaction.getId(), ga))
            return None
    return ''


def _remove_note_containing_text_of_interest(node, text):
    if not node:
        return False
    to_remove = []
    for i in range(0, node.getNumChildren()):
        child = node.getChild(i)
        note = child.getCharacters()
        if note and text in note \
                or _remove_note_containing_text_of_interest(child, text) and 0 == child.getNumChildren():
            to_remove.append(i)

    for i in reversed(to_remove):
        node.removeChild(i)
    return len(to_remove) > 0


def _get_nodes_of_type(node, type):
    if node:
        if node.isStart() and type in node.getName():
            yield node
        for i in range(0, node.getNumChildren()):
            for res in _get_nodes_of_type(node.getChild(i), type):
                yield res


def set_gene_association(reaction, gene_association):
    """Sets the reaction gene association. (The old gene association will be overwritten)
    """
    remove_gene_association(reaction)

    if not gene_association:
        return
    # Do not forget to convert s to str as if it is for example in Unicode, libsbml will complain
    s = str('<p>GENE_ASSOCIATION: %s</p>' % gene_association)

    body = next(_get_nodes_of_type(reaction.getNotes(), 'body'), None)
    if body:
        body.addChild(libsbml.XMLNode_convertStringToXMLNode(s))
    else:
        reaction.setNotes(libsbml.XMLNode_convertStringToXMLNode(
            str("<body xmlns='http://www.w3.org/1999/xhtml'>%s</body>" % s)))


def remove_gene_association(reaction):
    """Removes the reaction gene association.
    """
    _remove_note_containing_text_of_interest(reaction.getNotes(), GA_PREFIX)


def get_genes(reaction):
    """
    Extracts a set of genes from the gene association encoded as a note prefixed with GENE_ASSOCIATION:.
    :param reaction: the reaction of interest (libsbml.Reaction)
    :return: a set of genes of interest.
    """
    genes = set()

    def gene_parse_action(t):
        genes.add(t[0])
        return t[0]

    get_gene_association(reaction, gene_parse_action=gene_parse_action, flatten=False)
    return genes


def get_pathway_expression(reaction):
    result = set()
    node = reaction.getNotes()
    _get_prefixed_notes_value(node, result, PATHWAY_PREFIX)
    _get_prefixed_notes_value(node, result, "Pathway:")
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
    return {model.getSpecies(s_id).getCompartment() for s_id in get_metabolites(r, include_modifiers=False)}


def create_species(model, compartment_id, name=None, bound=False, id_=None, type_id=None, sbo_id=None):
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
    return new_species


def create_compartment(model, name=None, outside=None, id_=None):
    new_comp = model.createCompartment()
    id_ = generate_unique_id(model, id_ if id_ else "c")
    res = new_comp.setId(id_)
    if libsbml.LIBSBML_OPERATION_SUCCESS != res:
        logging.error("compartment %s creation error: %d" % (id_, res))
    if name:
        new_comp.setName(name)
    if outside:
        new_comp.setOutside(outside)
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
        if not id_[0].isalpha():
            id_ = 's_' + id_
        id_ = id_.encode('ascii', errors='ignore').decode()
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

    # for s in sorted((model.getSpecies(s_id) for s_id in s_ids_processed), key=lambda s: s.name):
    #     print(s.name, s.id, len(s_id2r_ids[s.id]))
    #
    # print("________________________________________")
    # for s in sorted((model.getSpecies(s_id) for s_id in ubiquitous_s_ids), key=lambda s: s.name):
    #     print(s.name, s.id)

    return r_ids


def create_reaction(model, r_id2st, p_id2st, name=None, reversible=True, id_=None):
    new_r_id = generate_unique_id(model, id_=id_)
    new_r = model.createReaction()
    if libsbml.LIBSBML_OPERATION_SUCCESS != new_r.setId(new_r_id):
        logging.error("reaction %s creation error" % new_r_id)
    if name:
        new_r.setName(name)
    new_r.setReversible(reversible)
    for m_id, st in r_id2st.items():
        sr = new_r.createReactant()
        sr.setSpecies(m_id)
        sr.setStoichiometry(st)
    for m_id, st in p_id2st.items():
        sr = new_r.createProduct()
        sr.setSpecies(m_id)
        sr.setStoichiometry(st)
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
