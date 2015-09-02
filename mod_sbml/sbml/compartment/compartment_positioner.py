from mod_sbml.annotation.rdf_annotation_helper import get_is_annotations, get_is_vo_annotations

GO_CYTOPLASM = 'go:0005737'
GO_CYTOSOL = 'go:0005829'
GO_CYTOPLASM_VARIANTS = {GO_CYTOPLASM, GO_CYTOSOL}
GO_NUCLEUS = 'go:0005634'
GO_ORGANELLE_OUTER_MEMBRANE = 'go:0031968'
GO_ORGANELLE_INNER_MEMBRANE = 'go:0019866'
GO_MEMBRANE = 'go:0016020'
GO_ENVELOPE = 'go:0031975'
GO_EXTRACELLULAR = 'go:0005576'
GO_ORGANELLE = 'go:0043226'
GO_CELL = 'go:0005623'
GO_LIPID_PARTICLE = 'go:0005811'

__author__ = 'anna'

partOfCheck = lambda it, t_id, onto: onto.part_of(it, [t_id])
partOf = lambda t_id, onto, candidate_ids: {it for it in candidate_ids if partOfCheck(it, t_id, onto)}
isACheck = lambda it, t_id, onto: onto.is_a(it, t_id) or (t_id.lower() == it.lower())
isA = lambda t_id, onto, candidate_ids: {it for it in candidate_ids if isACheck(it, t_id, onto)}
isAorPartOf = lambda t_id, onto, candidate_ids: {it for it in candidate_ids if
                                                 isACheck(it, t_id, onto) or partOfCheck(it, t_id, onto)}


def get_go_term(comp, onto):
    """
    Find a Gene Ontology term for a compartment of interest,
    using its annotations or name.
    :param comp: libSBML Compartment
    :param onto: the Gene Ontology
    :return: term if it was found otherwise None
    """
    for is_annotation in get_is_annotations(comp):
        term = onto.get_term(is_annotation, check_only_ids=False)
        if term:
            return term
    for is_vo_annotation in get_is_vo_annotations(comp):
        term = onto.get_term(is_vo_annotation, check_only_ids=False)
        if term:
            return term
    return onto.get_term(comp.getName() if comp.getName() else comp.getId(), check_only_ids=False)


def get_comp2go(model, onto):
    """
    Map compartments in the model to the Gene Ontology terms,
    using annotations in the model and compartment names.
    :param model: libSBML Model
    :param onto: the Gene Ontology
    :return: dict: comp_id -> go_term_id
    """
    comp2go = {}
    for comp in model.getListOfCompartments():
        term = get_go_term(comp, onto)
        if term:
            comp2go[comp.getId()] = term.get_id()
    return comp2go


def comp2level(model, onto):
    """
    Calculates levels of compartments in a given model, where 0 level corresponds to the outmost compartment(s),
    level i to the compartments surrounded by compartments of level i-1.
    :param model: libSBML model containing the compartments to be classified
    :param onto: the Gene Ontology
    :return: dict: compartment_id -> level, e.g. {Boundary: 0, Cytosol: 1, Mitochondrion: 2, Peroxisome: 2}
    """
    outs_are_set = next((True for comp in model.getListOfCompartments() if comp.getOutside() is not None
                         and model.getCompartment(comp.getOutside())), False)

    # if the compartment hierarchy is not provided in the model,\
    # let's infer it based on the Gene Ontology
    if not outs_are_set:
        term_id2comp_id = {}
        for comp in model.getListOfCompartments():
            term = get_go_term(comp, onto)
            if term:
                term_id2comp_id[term.get_id()] = comp.getId()
        in2out = nest_compartments_with_gene_ontology(set(term_id2comp_id.iterkeys()), onto)
        for in_t_id, out_t_id in in2out.iteritems():
            if out_t_id:
                model.getCompartment(term_id2comp_id[in_t_id]).setOutside(term_id2comp_id[out_t_id])

    roots = {comp.getId() for comp in model.getListOfCompartments() if not comp.isSetOutside()}
    outsides = {comp.getOutside() for comp in model.getListOfCompartments() if comp.isSetOutside()}

    # If there are several roots and one of them does not have any compartment inside it,
    # let's make it the global root.
    # (It is often a Boundary compartment that is not found in GO.)
    if len(roots) > 1:
        the_root = None
        for it in roots:
            if not (it in outsides):
                the_root = it
                break
        if the_root:
            for root in roots:
                if not root == the_root:
                    model.getCompartment(root).setOutside(the_root)
            roots = {the_root}

    i = 0
    c_id2level = {}
    current = roots
    while current:
        c_id2level.update({c_id: (i, model.getCompartment(c_id).getOutside()) for c_id in current})
        i += 1
        current = {comp.getId() for comp in model.getListOfCompartments() if
                   (comp.isSetOutside() and comp.getOutside() in current)}

    return c_id2level


def classify_cellular_components(t_ids, onto):
    """
    Classifies cellular components using the Gene Ontology
    :param t_ids: the Gene Ontology term ids to be classified
    :param onto: the Gene Ontology
    :return: dict that maps compartment_id to outside_compartment_id
    """

    # cell
    inside_cell = isAorPartOf(GO_CELL, onto, t_ids)

    # organelles and their parts
    organelle_parts = partOf(GO_ORGANELLE, onto, t_ids) | partOf(GO_NUCLEUS, onto, t_ids) \
                      | partOf(GO_LIPID_PARTICLE, onto, t_ids)

    not_org_parts = t_ids - organelle_parts
    organelles = isA(GO_ORGANELLE, onto, not_org_parts) | isA(GO_NUCLEUS, onto, not_org_parts) | \
                 isA(GO_LIPID_PARTICLE, onto, not_org_parts)
    organelle_parts |= organelles

    inside_cell |= organelle_parts
    outside_cell = t_ids - inside_cell

    organelle2parts = {it: set() for it in organelles}
    organelle_parts_wo_organelle = []
    for it in organelle_parts:
        orgs = onto.part_of(it, organelles)
        if orgs:
            organelle2parts[orgs.pop()].add(it)
        elif it not in organelle2parts:
            organelle_parts_wo_organelle.append(it)

    # find organelles for which only parts are present in our model
    organelle_ids = onto.get_descendants(GO_ORGANELLE, False) | {GO_ORGANELLE.lower()}
    for it in organelle_parts_wo_organelle:
        organelle = get_outside_comp_id(it, onto, organelle_ids)
        if organelle:
            if organelle in organelle2parts:
                organelle2parts[organelle].add(it)
            else:
                # (*) add an additional organelle
                organelle2parts[organelle] = {it}

    # extracellular
    extracellular = isAorPartOf(GO_EXTRACELLULAR, onto, outside_cell)

    return organelle2parts, extracellular, inside_cell, outside_cell


def classify_organelle_parts(parts, onto):
    """
    Classifies organelle parts into outer membranes, other membranes, other envelopes, inner membranes, organelles,
    and parts located inside organelles
    :param parts: ids of the GO terms of organelle parts to be classified
    :param onto: the Gene Ontology (GO)
    :return: tuple of sets of GO ids corresponding to (outer membranes, other membranes, other envelopes,
    inner membranes, organelles, parts located inside organelles)
    """
    # envelope, membrane
    envelopes = isAorPartOf(GO_ENVELOPE, onto, parts) | isAorPartOf(GO_MEMBRANE, onto, parts)
    insides = parts - envelopes
    # organelles
    organelles = isA(GO_ORGANELLE, onto, insides)
    insides -= organelles
    # membranes
    other_membranes = isA(GO_MEMBRANE, onto, envelopes)
    other_envelopes = envelopes - other_membranes
    # organelle inner membranes
    in_membranes = isA(GO_ORGANELLE_INNER_MEMBRANE, onto, other_membranes)
    other_membranes -= in_membranes
    # organelle outer membranes
    out_membranes = isA(GO_ORGANELLE_OUTER_MEMBRANE, onto, other_membranes)
    other_membranes -= out_membranes

    return out_membranes, other_membranes, other_envelopes, in_membranes, organelles, insides


def nest_compartments_with_gene_ontology(t_ids, onto):
    """
    Update compartment hierarchy using the Gene Ontology
    :param t_ids: the Gene Ontology term ids for which the hierarchy is to be found
    :param onto: the Gene Ontology
    :return: dict that maps compartment_term_id to outside_compartment_term_id
    """
    # Look for outside compartments in the Gene Ontology
    comp2out = {t_id: get_outside_comp_id(t_id, onto, t_ids) for t_id in t_ids}

    organelle2parts, extracellular, inside_cell, outside_cell = classify_cellular_components(t_ids, onto)
    # we might have added organelles that were not among t_ids, so let's add outsides for them
    comp2out.update({organelle: get_outside_comp_id(organelle, onto, t_ids) for organelle in organelle2parts.iterkeys()})

    organelle_parts = reduce(lambda s1, s2: s1 | s2, organelle2parts.itervalues(), set(organelle2parts.iterkeys()))

    innermost_extracellular_comp = None
    if extracellular:
        innermost_extracellular_comp, _ = \
            correct_membranes(None, extracellular, comp2out, onto)

    outermost_cell_part = get_outer_most(onto, inside_cell - organelle_parts)
    cell_innermost, cell_outermost = \
        correct_membranes(outermost_cell_part, (inside_cell - organelle_parts) - {outermost_cell_part},
                          comp2out, onto)

    # surround those that are not surrounded
    for comp in comp2out.iterkeys():
        if comp in inside_cell:
            if not comp2out[comp]:
                if comp in organelle_parts and comp != cell_outermost:
                    comp2out[comp] = cell_outermost
                else:
                    comp2out[comp] = innermost_extracellular_comp

    # correct organelle_membrane part_of organelle inferences:
    for organelle, parts in organelle2parts.iteritems():
        correct_membranes(organelle, parts, comp2out, onto)

    # remove additional terms added before in (*)
    additional_organelles = (set(comp2out.iterkeys()) | {it for it in comp2out.itervalues() if it}) - t_ids
    for it in additional_organelles:
        outside = comp2out[it]
        comp2out = {comp: (out if it != out else outside) for (comp, out) in comp2out.iteritems() if comp != it}

    return comp2out


def correct_membranes(organelle, parts, comp2out, onto):
    """
    This method ensures the following element inclusion, using the Gene Ontology:
    the compartment that was specified as outside of the organelle (if any) in comp2out is_outside_of
       organelle_outer_membranes are_outside_of
           other_organelle_membranes are_outside_of
               other_envelopes are_outside_of
                   organelle_inner_membranes are_outside_of
                       organelles are_outside_of
                           inside_parts_of_organelles
    :param organelle: the GO term of an organelle of interest (if any)
    :param parts: the GO term of a organelle parts of interest (if any)
    :param comp2out: current inclusion hierarchy, that will be updated in this method
    :param onto: the Gene Ontology
    :return: the innermost element and the outermost element
    """

    def inside(ins, outs):
        while True:
            oo = comp2out[ins]
            if not oo or oo == ins:
                return False
            if oo == outs:
                return True
            ins = oo

    # process envelopes
    def process(items, cur_outer, outer, inner):
        if items:
            if cur_outer:
                comp2out.update({it: cur_outer for it in items})
            cur_outer = next(iter(items))
            if not outer:
                outer = cur_outer
            inner = cur_outer
        return cur_outer, inner, outer

    outside_comp = comp2out[organelle] if organelle \
        else next((comp2out[it] for it in parts if not comp2out[it] in parts), None)

    current_outer = outside_comp
    innermost, outermost = None, None

    out_membranes, other_membranes, other_envelopes, in_membranes, organelles, insides = \
        classify_organelle_parts(parts, onto)

    for elements in (out_membranes, other_membranes, other_envelopes, in_membranes):
        current_outer, innermost, outermost = \
            process(elements, current_outer, innermost, outermost)

    # b) process organelles
    if current_outer:
        if organelle:
            comp2out[organelle] = current_outer
        for it in organelles:
            comp2out[it] = current_outer
    if organelle:
        current_outer = organelle
        if not outermost:
            outermost = current_outer
        innermost = current_outer

    # c) process insides
    if current_outer:
        for it in insides:
            if comp2out[it] not in insides:
                comp2out[it] = current_outer
            elif not inside(innermost, it):
                innermost = it
    else:
        for it in insides:
            if not innermost:
                innermost = it
            if not outermost:
                outermost = it
            if comp2out[it] not in insides:
                outermost = it
            elif not inside(innermost, it):
                innermost = it
    return innermost, outermost


def get_outside_comp_id(comp_id, onto, options):
    """
    Find the compartment among the given options (if it exists)
    that is the closest from outside to the compartment of interest.
    :param comp_id: the Gene Ontology term id of the compartment of interest.
    :param onto: the Gene Ontology
    :param options: the Gene Ontology term ids of the potential outside compartments
    :return: the Gene Ontology term id of the compartment among the given options
    that is the closest from outside to the compartment of interest, if it exists, None otherwise.
    """
    options = {it.lower() for it in options}
    # we will surround nucleus and lipid particle by cytoplasm manually,
    # as it is not a membrane-bounded organelle
    # => not part of the cytoplasm according to the Gene Ontology
    if isACheck(comp_id, GO_NUCLEUS, onto) or isACheck(comp_id, GO_LIPID_PARTICLE, onto):
        cyto_intersection = GO_CYTOPLASM_VARIANTS & options
        if cyto_intersection:
            return next(iter(cyto_intersection))
    candidates = set(options)
    candidates -= {comp_id}
    matches = onto.part_of(comp_id, candidates)
    return get_inner_most(onto, matches)


def get_inner_most(onto, options):
    """
    Find the innermost compartment among the given options.
    :param onto: the Gene Ontology
    :param options: the Gene Ontology term ids of the compartments of interest
    :return: the Gene Ontology term id of the innermost compartment among the given options.
    """
    if not options:
        return None
    if len(options) == 1:
        return options.pop()
    # return id of the inner-most compartment
    while options:
        it = options.pop()
        no_better_candidate = True
        for m in options:
            if onto.is_a(m, it) or onto.part_of(m, [it]):
                no_better_candidate = False
                break
        if no_better_candidate:
            return it
    return None


def get_outer_most(onto, options):
    """
    Find the outermost compartment among the given options.
    :param onto: the Gene Ontology
    :param options: the Gene Ontology term ids of the compartments of interest
    :return: the Gene Ontology term id of the outermost compartment among the given options.
    """
    if not options:
        return None
    if len(options) == 1:
        return options.pop()
    # return id of the outer-most compartment
    while options:
        it = options.pop()
        no_better_candidate = True
        for m in options:
            if onto.is_a(it, m) or onto.part_of(it, [m]):
                no_better_candidate = False
                break
        if no_better_candidate:
            return it
    return None

