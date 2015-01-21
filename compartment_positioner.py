from libsbml import BQB_IS, BQB_IS_VERSION_OF

from obo_ontology import miriam_to_term_id
from sbml_manager import get_qualifier_values, create_compartment


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
partOf = lambda t_id, onto, candidates: {it for it in candidates if partOfCheck(it.get_id(), t_id, onto)}
isACheck = lambda it, t_id, onto: onto.is_a(it, onto.get_term(t_id)) or (t_id.lower() == it.get_id().lower())
isA = lambda t_id, onto, candidates: {it for it in candidates if isACheck(it, t_id, onto)}
isAorPartOf = lambda t_id, onto, candidates: {it for it in candidates if
                                              isACheck(it, t_id, onto) or partOfCheck(it.get_id(), t_id, onto)}


def get_go_term(annotation, qualifier, onto):
    if not annotation:
        return None
    for go_id in get_qualifier_values(annotation, qualifier):
        go_id = miriam_to_term_id(go_id)
        term = onto.get_term(go_id)
        if term:
            return term
    return None


def get_comp2go(model, onto):
    comp2go = {}
    for comp in model.getListOfCompartments():
        annotation = comp.getAnnotation()
        term = get_go_term(annotation, BQB_IS, onto)
        if not term:
            term = get_go_term(annotation, BQB_IS_VERSION_OF, onto)
        if not term:
            term_ids = onto.get_ids_by_name(comp.getName())
            if term_ids:
                term = onto.get_term(term_ids.pop())
            else:
                term = None
        comp2go[comp.getId()] = term.get_id() if term else ''
    return comp2go


def comp2level(model, onto):
    outs_set = False
    for comp in model.getListOfCompartments():
        out = comp.getOutside()
        if out:
            out = model.getCompartment(out)
            if out:
                outs_set = True
    if not outs_set:
        term2comp = {}
        for comp in model.getListOfCompartments():
            annotation = comp.getAnnotation()
            if annotation:
                term = get_go_term(annotation, BQB_IS, onto)
                if not term:
                    term = get_go_term(annotation, BQB_IS_VERSION_OF, onto)
                if term:
                    term2comp[term] = comp
                    continue
            term_ids = onto.get_ids_by_name(comp.getName())
            if term_ids:
                term2comp[onto.get_term(set(term_ids).pop())] = comp
        in2out = nest_compartments_with_gene_ontology({it.get_id() for it in term2comp.iterkeys()}, onto)
        for in_term, out_term in in2out.iteritems():
            if out_term:
                term2comp[in_term].setOutside(term2comp[out_term].getId())
    roots = {comp.getId() for comp in model.getListOfCompartments() if not comp.isSetOutside()}
    outsides = {comp.getOutside() for comp in model.getListOfCompartments() if comp.isSetOutside()}
    if len(roots) > 1:
        the_root = None
        for it in roots:
            if not (it in outsides):
                the_root = it
                break
        if not the_root:
            comp_name = model.getName()
            if not comp_name:
                comp_name = model.getId()
            the_root = create_compartment(model, comp_name).getId()
        for root in roots:
            if not root == the_root:
                model.getCompartment(root).setOutside(the_root)
    else:
        the_root = roots.pop()
    i = 0
    c_id2level = {}
    current = {the_root}
    while current:
        c_id2level.update({c_id: (i, model.getCompartment(c_id).getOutside()) for c_id in current})
        i += 1
        current = {comp.getId() for comp in model.getListOfCompartments() if
                   (comp.isSetOutside() and comp.getOutside() in current)}

    return c_id2level


# Update compartment hierarchy using the Gene Ontology
def nest_compartments_with_gene_ontology(t_ids, onto):
    comp2out = {}
    # Look for missing outside compartments in the Gene Ontology
    for t_id in t_ids:
        term = onto.get_term(t_id)
        if not term:
            continue
        out = get_outside_comp_id(t_id, onto, t_ids)
        if out:
            out = onto.get_term(out)
        comp2out[term] = out

    comps = set(comp2out.iterkeys())
    # cell
    inside_cell = isAorPartOf(GO_CELL, onto, comps)

    # organelle
    organelle_parts = partOf(GO_ORGANELLE, onto, comps) | partOf(GO_NUCLEUS, onto, comps) | partOf(GO_LIPID_PARTICLE,
                                                                                                   onto, comps)

    not_org_parts = comps - organelle_parts
    organelles = isA(GO_ORGANELLE, onto, not_org_parts) | isA(GO_NUCLEUS, onto, not_org_parts) | isA(GO_LIPID_PARTICLE,
                                                                                                     onto,
                                                                                                     not_org_parts)
    organelle_parts |= organelles

    inside_cell |= organelle_parts
    outside_cell = comps - inside_cell

    # extracellular
    extracellulars = isAorPartOf(GO_EXTRACELLULAR, onto, outside_cell)
    extracellular_ids = {c.get_id() for c in extracellulars}
    extracellular_comp = None
    if extracellular_ids:
        extracellular_comp, eo = correct_membranes(None, extracellulars, comp2out, onto)

    outer = onto.get_term(get_outer_most(onto, [it.get_id() for it in inside_cell - organelle_parts]))
    cell_inner, cell_outer = correct_membranes(outer, (inside_cell - organelle_parts) - {outer}, comp2out, onto)

    organelle2parts = {it: set() for it in organelles}
    organelle_ids = {it.get_id() for it in organelles}
    no_organelle_parts = []
    for it in organelle_parts:
        orgs = onto.part_of(it.get_id(), organelle_ids)
        if orgs:
            organelle2parts[onto.get_term(orgs.pop())].add(it)
        else:
            no_organelle_parts.append(it)

    populated = {org for org in organelle2parts.iterkeys() if lambda org: organelle2parts[org]}
    # organelle
    organelle = onto.get_term(GO_ORGANELLE)
    organelle_ids = {it.get_id() for it in organelle.get_descendants(False)} | {GO_ORGANELLE.lower()}
    for it in no_organelle_parts:
        if it in organelles:
            if it in populated:
                continue
            parents = set(onto.get_descendants(it, False)) & populated
            if parents:
                organelle2parts[parents.pop()].add(it)
                continue
        organelle = get_outside_comp_id(it.get_id(), onto, organelle_ids)
        if organelle:
            if organelle in organelle2parts:
                organelle2parts[organelle].add(it)
            else:
                organelle2parts[organelle] = {it}
                out = get_outside_comp_id(organelle, onto, t_ids)
                if out:
                    out = onto.get_term(out)
                comp2out[organelle] = out

    # surround those that are not surrounded, by extracellular
    for comp in comp2out.iterkeys():
        if comp in inside_cell:
            if not comp2out[comp]:
                if comp in organelle_parts and comp != cell_outer:
                    comp2out[comp] = cell_outer
                else:
                    comp2out[comp] = extracellular_comp

    # Correct organelle_membrane part_of organelle inferences:
    for organelle, parts in organelle2parts.iteritems():
        correct_membranes(organelle, parts, comp2out, onto)
    return comp2out


def correct_membranes(organelle, parts, comp2out, onto):
    outside = None
    if organelle:
        outside = comp2out[organelle]
    else:
        for it in parts:
            out = comp2out[it]
            if not (out in parts):
                outside = out
                break
    # envelope, membrane
    envelopes = isAorPartOf(GO_ENVELOPE, onto, parts) | isAorPartOf(GO_MEMBRANE, onto, parts)
    insides = parts - envelopes
    # organelle
    organelles = isA(GO_ORGANELLE, onto, insides)
    insides -= organelles
    # membrane
    membranes = isA(GO_MEMBRANE, onto, envelopes)
    envelope_others = envelopes - membranes
    # organelle inner membrane
    in_membranes = isA(GO_ORGANELLE_INNER_MEMBRANE, onto, membranes)
    membranes -= in_membranes
    # organelle outer membrane
    out_membranes = isA(GO_ORGANELLE_OUTER_MEMBRANE, onto, membranes)
    membranes -= out_membranes

    inner, outer = None, None
    for it in envelopes:
        comp2out[it] = outside
    out = None
    if out_membranes:
        out = out_membranes.pop()
    if out:
        if not outer:
            outer = out
        inner = out
        for it in membranes:
            comp2out[it] = out
    if membranes:
        out = membranes.pop()
    if out:
        if not outer:
            outer = out
        inner = out
        for it in envelope_others:
            comp2out[it] = out
    if envelope_others:
        out = envelope_others.pop()
    if out:
        if not outer:
            outer = out
        inner = out
        for it in in_membranes:
            comp2out[it] = out
    if in_membranes:
        out = in_membranes.pop()
    if out:
        if not outer:
            outer = out
        inner = out
        if organelle:
            comp2out[organelle] = out
        for it in organelles:
            comp2out[it] = out
    if not isinstance(organelle, str):
        out = organelle
    else:
        del comp2out[organelle]

    def inside(ins, outs):
        while True:
            oo = comp2out[ins]
            if oo == ins:
                return False
            if not oo:
                return False
            if oo == outs:
                return True
            ins = oo

    if out:
        if not outer:
            outer = out
        inner = out
        for it in insides:
            if not (comp2out[it] in insides):
                comp2out[it] = out
            elif not inside(inner, it):
                inner = it
    else:
        for it in insides:
            if not inner:
                inner = it
            if not outer:
                outer = it
            if not (comp2out[it] in insides):
                outer = it
            elif not inside(inner, it):
                inner = it
    return inner, outer


def get_outside_comp_id(comp_id, onto, variants):
    # we will surround nucleus by cytoplasm manually,
    # as it is not a membrane-bounded organelle
    # => not part of the cytoplasm according to the Gene Ontology
    # nucleus, cytoplasm
    variants = {it.lower() for it in variants}
    if isACheck(onto.get_term(comp_id), GO_NUCLEUS, onto) and (GO_CYTOPLASM in variants):
        return GO_CYTOPLASM
    if isACheck(onto.get_term(comp_id), GO_LIPID_PARTICLE, onto) and (GO_CYTOPLASM in variants):
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
            if onto.is_a(onto.get_term(m), onto.get_term(it)) or onto.part_of(m, [it]):
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
            if onto.is_a(onto.get_term(it), onto.get_term(m)) or onto.part_of(it, [m]):
                no_better_candidate = False
                break
        if no_better_candidate:
            return it
    return None

