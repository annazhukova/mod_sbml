from mod_sbml.annotation.rdf_annotation_helper import get_is_annotations, get_is_vo_annotations

GO_CYTOPLASM = 'go:0005737'
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
    comp2go = {}
    for comp in model.getListOfCompartments():
        term = get_go_term(comp, onto)
        if term:
            comp2go[comp.getId()] = term.get_id()
    return comp2go


def comp2level(model, onto):
    outs_set = next((True for comp in model.getListOfCompartments() if comp.getOutside() is not None
                     and model.getCompartment(comp.getOutside())), False)
    if not outs_set:
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
    no_organelle_parts = []
    for it in organelle_parts:
        orgs = onto.part_of(it, organelles)
        if orgs:
            organelle2parts[orgs.pop()].add(it)
        elif it not in organelle2parts:
            no_organelle_parts.append(it)

    return organelle2parts, no_organelle_parts, inside_cell, outside_cell


def nest_compartments_with_gene_ontology(t_ids, onto):
    """
    Update compartment hierarchy using the Gene Ontology
    :param t_ids: the Gene Ontology term ids for which the hierarchy is to be found
    :param onto: the Gene Ontology
    :return: dict that maps compartment_term_id to outside_compartment_term_id
    """
    # Look for outside compartments in the Gene Ontology
    comp2out = {t_id: get_outside_comp_id(t_id, onto, t_ids) for t_id in t_ids}

    organelle2parts, no_organelle_parts, inside_cell, outside_cell = classify_cellular_components(t_ids, onto)
    organelle_parts = reduce(lambda s1, s2: s1 | s2, organelle2parts.itervalues(), set(organelle2parts.iterkeys()))

    # extracellular
    extracellulars = isAorPartOf(GO_EXTRACELLULAR, onto, outside_cell)
    innermost_extracellular_comp = None
    if extracellulars:
        innermost_extracellular_comp, _ = \
            correct_membranes(None, extracellulars, comp2out, onto)

    outermost_cell_part = get_outer_most(onto, inside_cell - organelle_parts)
    cell_innermost, cell_outermost = \
        correct_membranes(outermost_cell_part, (inside_cell - organelle_parts) - {outermost_cell_part},
                          comp2out, onto)

    # find organelles for which only parts are present in our model
    organelle_ids = onto.get_descendants(GO_ORGANELLE, False) | {GO_ORGANELLE.lower()}
    for it in no_organelle_parts:
        organelle = get_outside_comp_id(it, onto, organelle_ids)
        if organelle:
            if organelle in organelle2parts:
                organelle2parts[organelle].add(it)
            else:
                # (*) add an additional organelle
                organelle2parts[organelle] = {it}
                out_id = get_outside_comp_id(organelle, onto, t_ids)
                if out_id:
                    comp2out[organelle] = out_id

    # surround those that are not surrounded, by extracellular
    for comp in comp2out.iterkeys():
        if comp in inside_cell:
            if not comp2out[comp]:
                if comp in organelle_parts and comp != cell_outermost:
                    comp2out[comp] = cell_outermost
                else:
                    comp2out[comp] = innermost_extracellular_comp

    # Correct organelle_membrane part_of organelle inferences:
    for organelle, parts in organelle2parts.iteritems():
        correct_membranes(organelle, parts, comp2out, onto)

    # Remove additional organelle added in (*)
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
    outside_comp = comp2out[organelle] if organelle \
        else next((comp2out[it] for it in parts if not comp2out[it] in parts), None)
    # envelope, membrane
    envelopes = isAorPartOf(GO_ENVELOPE, onto, parts) | isAorPartOf(GO_MEMBRANE, onto, parts)
    insides = parts - envelopes
    # organelle
    organelles = isA(GO_ORGANELLE, onto, insides)
    insides -= organelles
    # membrane
    other_membranes = isA(GO_MEMBRANE, onto, envelopes)
    other_envelopes = envelopes - other_membranes
    # organelle inner membrane
    in_membranes = isA(GO_ORGANELLE_INNER_MEMBRANE, onto, other_membranes)
    other_membranes -= in_membranes
    # organelle outer membrane
    out_membranes = isA(GO_ORGANELLE_OUTER_MEMBRANE, onto, other_membranes)
    other_membranes -= out_membranes

    def inside(ins, outs):
        while True:
            oo = comp2out[ins]
            if not oo or oo == ins:
                return False
            if oo == outs:
                return True
            ins = oo

    current_outer = outside_comp
    innermost, outermost = None, None

    # a) process membranes
    def process(elements, cur_outer, outer, inner):
        if elements:
            if cur_outer:
                for it in elements:
                    comp2out[it] = cur_outer
            cur_outer = next(iter(elements))
            if not outer:
                outer = cur_outer
            inner = cur_outer
        return cur_outer, inner, outer

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


def get_outside_comp_id(comp_id, onto, variants):
    variants = {it.lower() for it in variants}
    # we will surround nucleus and lipid particle by cytoplasm manually,
    # as it is not a membrane-bounded organelle
    # => not part of the cytoplasm according to the Gene Ontology
    if isACheck(comp_id, GO_NUCLEUS, onto) and (GO_CYTOPLASM in variants):
        return GO_CYTOPLASM
    if isACheck(comp_id, GO_LIPID_PARTICLE, onto) and (GO_CYTOPLASM in variants):
        return GO_CYTOPLASM
    candidates = set(variants)
    candidates -= {comp_id}
    matches = onto.part_of(comp_id, candidates)
    return get_inner_most(onto, matches)


def get_inner_most(onto, matches):
    if not matches:
        return None
    if len(matches) == 1:
        return matches.pop()
    # return id of the inner-most compartment
    while matches:
        it = matches.pop()
        no_better_candidate = True
        for m in matches:
            if onto.is_a(m, it) or onto.part_of(m, [it]):
                no_better_candidate = False
                break
        if no_better_candidate:
            return it
    return None


def get_outer_most(onto, matches):
    if not matches:
        return None
    if len(matches) == 1:
        return matches.pop()
    # return id of the outer-most compartment
    while matches:
        it = matches.pop()
        no_better_candidate = True
        for m in matches:
            if onto.is_a(it, m) or onto.part_of(it, [m]):
                no_better_candidate = False
                break
        if no_better_candidate:
            return it
    return None

