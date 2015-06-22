from collections import defaultdict
import re

import libsbml
import openpyxl

from openpyxl.styles import Style, Font

from mod_sbml.gibbs.gibbs_manager import get_r_kegg2equilibrium_const
from mod_sbml.gibbs.reaction_boundary_manager import get_bounds, compare_kegg_formulas, NO_KEGG_ID, NO_KEGG_INFO_FOUND, \
    get_better_r_bounds, REVERSED, SUCCESS
from mod_sbml.kegg.reaction_manager import get_rn2kegg_formula, get_rn2compounds
from mod_sbml.sbml.sbml_manager import get_kegg_m_id, get_kegg_r_id, get_gene_association, get_r_comps, get_reactants, \
    get_products, get_metabolites
from mod_sbml.serialization.xlsx_helper import save_data, add_values, HEADER_STYLE, RED_STYLE, BASIC_STYLE

__author__ = 'anna'


def serialize_model_info(sbml, path, r_style=lambda r_id: BASIC_STYLE):
    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()
    wb = openpyxl.Workbook()
    save_data(["Id", "Name", "Compartment", "Kegg"], (
        (m.id, m.name, model.getCompartment(m.getCompartment()).name, get_kegg_m_id(m))
        for m in model.getListOfSpecies()), ws_index=0, ws_name="Metabolites", wb=wb)
    save_data(["Id", "Name", "Lower Bound", "Upper Bound", "Formula", "Kegg", "Gene association"],
              ((r.id, r.name, get_bounds(r)[0], get_bounds(r)[1], get_sbml_r_formula(model, r, False), get_kegg_r_id(r),
                get_gene_association(r))
               for r in model.getListOfReactions()), styles=(r_style(r.id) for r in model.getListOfReactions()),
              ws_index=1, ws_name="Reactions", wb=wb)
    wb.save(path)


def serialize_info(model, r_ids, path):
    m_ids = set()
    wb = openpyxl.Workbook()
    ws = wb.create_sheet(1, "Reactions")
    add_values(ws, 1, 1,
               ["Id", "Name", "Low. Bound", "Upp. Bound", "Formula", "KEGG Formula", "Kegg", "Genes",
                "Ke(pH=7,T=298.15K)"],
               HEADER_STYLE)
    row = 2
    kegg2formula = get_rn2kegg_formula()
    k_formula = lambda kegg: kegg2formula[kegg] if kegg and kegg in kegg2formula else ''
    kegg2ke = get_r_kegg2equilibrium_const(ph=7)
    ke = lambda kegg: kegg2ke[kegg] if kegg in kegg2ke else ''
    for r_id in r_ids:
        r = model.getReaction(r_id)
        lb, ub = get_bounds(r)
        kegg = get_kegg_r_id(r)
        values = [r.id, r.name, lb, ub, get_sbml_r_formula(model, r, False), k_formula(kegg), kegg,
                  get_gene_association(r), ke(kegg)]
        add_values(ws, row, 1, values)
        row += 1
        m_ids |= get_metabolites(r)
    save_data(["Id", "Name", "Compartment", "Kegg"],
              ((m.id, m.name, model.getCompartment(m.getCompartment()).name, get_kegg_m_id(m))
               for m in (model.getSpecies(m_id) for m_id in m_ids)),
              ws_index=0, ws_name="Metabolites", wb=wb)
    wb.save(path)


def serialize_info_with_pw(model, pw2r_ids, o_r_ids, pw2pw_name, path):
    wb = openpyxl.Workbook()
    save_data(["Id", "Name", "Compartment", "Kegg"],
              ((m.id, m.name, model.getCompartment(m.getCompartment()).name, get_kegg_m_id(m))
               for m in model.getListOfSpecies()),
              ws_index=0, ws_name="Metabolites", wb=wb)
    if o_r_ids:
        pw2r_ids["No pathway"] = o_r_ids
    kegg2ke = get_r_kegg2equilibrium_const(ph=7)
    kegg2formula = get_rn2kegg_formula()
    kegg2rs_ps = get_rn2compounds()
    m2kegg = {m.id: get_kegg_m_id(m) for m in model.getListOfSpecies()}
    k_formula = lambda kegg: kegg2formula[kegg] if kegg and kegg in kegg2formula else ''
    for pw in sorted(pw2r_ids.iterkeys(), key=pw2pw_name):
        r_ids = pw2r_ids[pw]
        if not r_ids:
            continue
        pw_name = pw2pw_name(pw)
        pw = pw_name if pw_name else pw
        pw = pw[:min(31, len(pw))]
        pw = ''.join(e for e in pw if e.isalnum() or e.isspace())
        # ws = wb.create_sheet(title=pw)
        # add_values(ws, 1, 1,
        # ["Id", "Name", "Formula", "Kegg Formula", "Genes", "Kegg", "Ke(pH=7,T=298.15K)", "Low. B.", "Upp. B.", "New L. B.", "New Up. B.", "Comment"],#, "KEGG Formula", "Ke(pH=7,T=298.15K)"],
        #            HEADER_STYLE)
        # row = 2
        data, styles = [], []
        for r_id in sorted(r_ids):
            r = model.getReaction(r_id)
            lb, ub = get_bounds(r)
            kegg = get_kegg_r_id(r)
            result = compare_kegg_formulas(r, kegg2rs_ps, m2kegg)
            r_lower_bound, r_upper_bound = lb, ub
            ke = None
            if result in [NO_KEGG_ID, NO_KEGG_INFO_FOUND]:
                comment = "No KEGG formula found"
            elif kegg in kegg2ke and result in [SUCCESS, REVERSED]:
                ke = kegg2ke[kegg]
                r_lower_bound, r_upper_bound = get_better_r_bounds(r, ke, result == REVERSED, 1000)
                if result == REVERSED:
                    comment = "Reversed KEGG formula"
                else:
                    comment = ""
            else:
                comment = "Different to KEGG formula"

            style = RED_STYLE if lb != r_lower_bound or ub != r_upper_bound else BASIC_STYLE
            styles.append(style)
            data.append(
                [r.id, r.name, get_sbml_r_formula(model, r, False), k_formula(kegg), get_gene_association(r), kegg, ke,
                 lb, ub, r_lower_bound, r_upper_bound, comment])
            # add_values(ws, row, 1, values)
            # row += 1
        save_data(["Id", "Name", "Formula", "Kegg Formula", "Genes", "Kegg", "Ke(pH=7,T=298.15K)", "Low. B.", "Upp. B.",
                   "New L. B.", "New Up. B.", "Comment"],
                  data=data, styles=styles, ws_name=pw, wb=wb)
    wb.save(path)


def read_reactions(path):
    workbook = openpyxl.load_workbook(path)
    worksheet = workbook.get_active_sheet()
    metabolites = defaultdict(set)
    compartments = set()
    r_name2info = {}
    for row in worksheet.rows:
        r_name = str(row[1].value).encode('utf-8')
        r_formula = str(row[2].value).encode('utf-8')
        reversible = str(row[3].value).encode('utf-8').lower().find("irreversible") == -1
        reactants, products = r_formula.split(" = ")
        reactants = [r.replace('"', '').strip() for r in reactants.split(" + ")] if reactants.strip() else []
        products = [r.replace('"', '').strip() for r in products.split(" + ")] if products.strip() else []
        r2stoich = get_stoichiometries(reactants, metabolites, compartments)
        p2stoich = get_stoichiometries(products, metabolites, compartments)
        r_name2info[r_name] = r2stoich, p2stoich, reversible
    return metabolites, compartments, r_name2info


def get_stoichiometries(participants, metabolites, compartments):
    m2stoich = {}
    for m in participants:
        if m.find(" * ") != -1:
            st, m = m.split(" * ")
            st = int(st)
        else:
            st = 1
        comp = None
        comps = re.findall(r' \[\w+\]', m)
        for comp in comps:
            m = m.replace(comp, "")
            comp = comp[2:-1].strip()
            compartments.add(comp)
            break
        metabolites[m].add(comp)
        m2stoich[m] = st, comp
    return m2stoich


def get_sbml_r_formula(model, r, comp=True, id=True):
    format_m = lambda m_id, st: "%s%s" % (
        "%4.2f " % st if st != 1 else "",
        format_m_name(model.getSpecies(m_id), model, comp, id))
    formula = " + ".join([format_m(m_id, st) for (m_id, st) in get_reactants(r, True)]) + \
              (" <=> " if r.getReversible() else "-->") + \
              " + ".join([format_m(m_id, st) for (m_id, st) in get_products(r, True)])
    return formula


def format_m_name(m, model, comp=True, id=True):
    name = m.getName()
    c = model.getCompartment(m.getCompartment())
    c_name = c.getName()
    if -1 == name.find(c_name) and comp:
        name = "%s[%s]" % (name, c_name)
    if not comp:
        name = name.replace("[%s]" % c_name, "").replace("[%s]" % c.getId(), "").strip()
    return "%s(%s)" % (name, m.id) if id else name


def get_cobra_r_formula(r, comp=True):
    if comp:
        return "{0} <=> {1}".format(
            " + ".join(
                ["%s%s[%s](%s)" % (('%g ' % -r.get_coefficient(m)) if r.get_coefficient(m) != -1 else '',
                                   m.name, m.compartment, m.id) for m in r.reactants]),
            " + ".join(["%s%s[%s](%s)" % (('%g ' % r.get_coefficient(m)) if r.get_coefficient(m) != 1 else '',
                                          m.name, m.compartment, m.id) for m in r.products]))
    else:
        return "{0} <=> {1}".format(
            " + ".join(
                ["%s%s(%s)" % (('%g ' % -r.get_coefficient(m)) if r.get_coefficient(m) != -1 else '',
                               m.name, m.id) for m in r.reactants]),
            " + ".join(["%s%s(%s)" % (('%g ' % r.get_coefficient(m)) if r.get_coefficient(m) != 1 else '',
                                      m.name, m.id) for m in r.products]))


def save_model_mappings(source_model, target_model,
                        m_id2m_id, m_id2m_id_diff_comp, s_m_ids_unmapped, m2kegg, m2chebi,
                        r_id2r_id, r_id2r_id_diff_comp, s_r_ids_unmapped, r2kegg, r2bounds,
                        filename):
    workbook = openpyxl.Workbook()
    m_ws = workbook.create_sheet(0, "Metabolite mapping")
    save_metabolite_mappings(m_ws, source_model, target_model,
                             m_id2m_id, m_id2m_id_diff_comp, s_m_ids_unmapped, m2kegg, m2chebi)
    r_ws = workbook.create_sheet(1, "Reaction mapping")
    save_reaction_mappings(r_ws, source_model, target_model,
                           r_id2r_id, r_id2r_id_diff_comp, s_r_ids_unmapped, r2kegg, r2bounds)
    workbook.save(filename)


def save_metabolite_mappings(ws, source_model, target_model, m_id2m_id, m_id2m_id_diff_comp,
                             s_m_ids_unmapped, m2kegg, m2chebi):
    headers = ["Id", "Name", "Compartment", "KEGG", "ChEBI"]
    add_values(ws, 1, 2, [source_model.name if source_model.name else source_model.id], HEADER_STYLE)
    add_values(ws, 2, 2, headers, HEADER_STYLE)
    add_values(ws, 1, 2 + len(headers), [target_model.name if target_model.name else target_model.id], HEADER_STYLE)
    add_values(ws, 2, 2 + len(headers), headers, HEADER_STYLE)

    ws.cell(row=3, column=1).value = "Matching"
    i = 4
    for s_m_id in sorted(m_id2m_id.iterkeys()):
        add_metabolite(ws, s_m_id, source_model, i, 2, m2kegg, m2chebi)
        t_m_id = m_id2m_id[s_m_id]
        add_metabolite(ws, t_m_id, target_model, i, 2 + len(headers), m2kegg, m2chebi)
        i += 1

    ws.cell(row=i, column=1).value = "Similar in different compartments"
    i += 1
    for s_m_id in sorted(m_id2m_id_diff_comp.iterkeys()):
        add_metabolite(ws, s_m_id, source_model, i, 2, m2kegg, m2chebi)
        t_m_id = m_id2m_id_diff_comp[s_m_id]
        add_metabolite(ws, t_m_id, target_model, i, 2 + len(headers), m2kegg, m2chebi)
        i += 1

    ws.cell(row=i, column=1).value = "Unmapped"
    i += 1
    for s_m_id in sorted(s_m_ids_unmapped):
        add_metabolite(ws, s_m_id, source_model, i, 2, m2kegg, m2chebi)
        i += 1


def save_reaction_mappings(ws, source_model, target_model, r_id2r_id, r_id2r_id_diff_comp,
                           s_r_ids_unmapped, r2kegg, r2bounds):
    headers = ["Id", "Name", "Compartments", "Formula", "KEGG", "Genes", "Lower bound", "Upper bound",
               "Ke (pH=7, T=298.15K)"]
    add_values(ws, 1, 2, [source_model.name if source_model.name else source_model.id], HEADER_STYLE)
    add_values(ws, 2, 2, headers, HEADER_STYLE)
    add_values(ws, 1, 2 + len(headers), [target_model.name if target_model.name else target_model.id], HEADER_STYLE)
    add_values(ws, 2, 2 + len(headers), headers, HEADER_STYLE)

    kegg2ke = get_r_kegg2equilibrium_const(ph=7)
    ke = lambda kegg: kegg2ke[kegg] if kegg in kegg2ke else ''

    ws.cell(row=3, column=1).value = "Matching"
    i = 4
    for s_r_id in sorted(r_id2r_id.iterkeys()):
        t_r_id = r_id2r_id[s_r_id]
        s = Style(font=Font(color=openpyxl.styles.colors.BLACK))
        if s_r_id in r2bounds and t_r_id in r2bounds and r2bounds[s_r_id] != r2bounds[t_r_id]:
            s = Style(font=Font(color=openpyxl.styles.colors.RED))
        add_reaction(ws, s_r_id, source_model, i, 2, r2kegg, r2bounds, ke, s)
        add_reaction(ws, t_r_id, target_model, i, 2 + len(headers), r2kegg, r2bounds, ke, s)
        i += 1

    ws.cell(row=i, column=1).value = "Similar in different compartments"
    i += 1
    for s_r_id in sorted(r_id2r_id_diff_comp.iterkeys()):
        t_r_id = r_id2r_id_diff_comp[s_r_id]
        s = Style(font=Font(color=openpyxl.styles.colors.BLACK))
        if s_r_id in r2bounds and t_r_id in r2bounds and r2bounds[s_r_id] != r2bounds[t_r_id]:
            s = Style(font=Font(color=openpyxl.styles.colors.RED))
        add_reaction(ws, s_r_id, source_model, i, 2, r2kegg, r2bounds, ke, s)
        add_reaction(ws, t_r_id, target_model, i, 2 + len(headers), r2kegg, r2bounds, ke, s)
        i += 1

    ws.cell(row=i, column=1).value = "Unmapped"
    i += 1
    for s_r_id in sorted(s_r_ids_unmapped):
        add_reaction(ws, s_r_id, source_model, i, 2, r2kegg, r2bounds, ke)
        i += 1


def add_metabolite(ws, m_id, model, row, col, m2kegg, m2chebi):
    m = model.getSpecies(m_id)
    c = model.getCompartment(m.getCompartment())
    values = [m.id, format_m_name(m, model, False), c.name, m2kegg(m), m2chebi(m_id, model)]
    add_values(ws, row, col, values)


def add_reaction(ws, r_id, model, row, col, r2kegg, r2bounds, ke,
                 style=Style(font=Font(color=openpyxl.styles.colors.BLACK))):
    r = model.getReaction(r_id)
    cs = ", ".join(sorted((model.getCompartment(c_id).name for c_id in get_r_comps(r_id, model))))
    kegg = r2kegg(r)
    values = [r.id, r.name, cs, get_sbml_r_formula(model, r), kegg, get_gene_association(r)]
    if r_id in r2bounds:
        u_b, l_b = r2bounds[r_id]
        values.append(l_b)
        values.append(u_b)
    ke = ke(kegg)
    if ke:
        values.append(ke)
    add_values(ws, row, col, values, style)


def add_genes_to_reaction_list(sbml, path, out_path):
    workbook = openpyxl.load_workbook(path)
    worksheet = workbook.get_active_sheet()

    worksheet.cell(row=1, column=7).value = "Gene association"
    worksheet.cell(row=1, column=8).value = "KEGG Reaction"
    input_doc = libsbml.SBMLReader().readSBML(sbml)
    model = input_doc.getModel()
    r_len = len(worksheet.rows)
    for i in xrange(2, r_len + 1):
        r_id = str(worksheet.cell(row=i, column=1).value)
        r = model.getReaction(r_id)
        worksheet.cell(row=i, column=7).value = get_gene_association(r)
        worksheet.cell(row=i, column=8).value = get_kegg_r_id(r)
    workbook.save(out_path)

