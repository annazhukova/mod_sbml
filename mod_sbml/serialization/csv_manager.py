import pandas
from pandas import DataFrame

from mod_sbml.annotation.chebi.chebi_annotator import get_chebi_id
from mod_sbml.annotation.gene_ontology.go_annotator import get_go_id
from mod_sbml.annotation.kegg.kegg_annotator import get_kegg_r_id, get_kegg_m_id
from mod_sbml.sbml.reaction_boundary_manager import get_bounds
from mod_sbml.sbml.sbml_manager import get_gene_association, get_formulas, get_pathway_expression, get_r_comps
from mod_sbml.serialization import get_sbml_r_formula, df2csv

__author__ = 'anna'


def serialize_model_info(model, prefix, c_id2level=None, blocked_reactions=[], to_excel=False):
    """
    Serializes the information on model compartments, metabolites and reactions into 3 tab-delimited files:
    <prefix>_compartments.tab, <prefix>_metabolites.tab, and <prefix>_reactions.tab.
    :param model: libsbml.Model the model to be serialized
    :param prefix: str, path prefix for the output files: <prefix>_compartments.tab, <prefix>_metabolites.tab,
    and <prefix>_reactions.tab
    :param c_id2level: (optional), ordering of compartment dict {c_id: level} into levels (int).
    If specified, reactions and metabolites will be sorted by compartment.
    :return: void
    """

    name2df = {}
    for (to_df, sheet) in ((compartments2df, 'compartments'), (metabolites2df, 'metabolites'),
                           (lambda model, c_id2level: reactions2df(model, c_id2level=c_id2level,
                                                                   blocked_reactions=blocked_reactions), 'reactions')):
        df = to_df(model, c_id2level=c_id2level)
        df.dropna(axis=1, how='all', inplace=True)
        name2df[sheet] = df

    if to_excel:
        with pandas.ExcelWriter('%smodel.xlsx' % prefix) as writer:
            for (sheet, df) in name2df.items():
                df.to_excel(writer, sheet_name=sheet, header=True, index=False, na_rep='')
    else:
        for (name, df) in name2df.items():
            csv = '%s%s.tab' % (prefix, name)
            df2csv(df, csv)


def metabolites2df(model, c_id2level=None):
    data = []
    index = []

    def get_key(m):
        c_id = m.getCompartment()
        if c_id2level:
            return c_id2level[c_id], c_id, m.name
        return m.name, c_id

    for m in sorted(model.getListOfSpecies(), key=get_key):
        formulas = get_formulas(m)
        data.append(
            (m.id, m.name, model.getCompartment(m.getCompartment()).getName(),
             formulas.pop() if formulas else None, get_kegg_m_id(m), get_chebi_id(m)))
        index.append(m.id)
    columns = ['Id', 'Name', 'Compartment', 'Formula', 'KEGG', 'ChEBI']
    return DataFrame(data=data, index=index, columns=columns)


def compartments2df(model, c_id2level=None):
    data = []
    index = []

    def get_key(c):
        if c_id2level:
            return c_id2level[c.id], c.name, c.id
        return c.name, c.id

    for c in sorted(model.getListOfCompartments(), key=get_key):
        data.append((c.id, c.name, get_go_id(c)))
        index.append(c.id)
    return DataFrame(data=data, index=index, columns=['Id', 'Name', "GO"])


def reactions2df(model, r_ids=None, c_id2level=None, blocked_reactions=None):
    if blocked_reactions is None:
        blocked_reactions = []
    data = []
    index = []

    def get_key(r):
        c_ids = tuple(sorted(get_r_comps(r.id, model)))

        if c_id2level:
            return tuple(sorted({c_id2level[c_id] for c_id in c_ids})), c_ids, r.id
        return r.id

    rs = model.getListOfReactions() if not r_ids else (r for r in model.getListOfReactions() if r.id in r_ids)

    for r in sorted(rs, key=get_key):
        record = (r.id, r.name, get_bounds(r)[0], get_bounds(r)[1],
                  get_sbml_r_formula(model, r, show_compartments=True, show_metabolite_ids=False),
                  ', '.join(tuple(sorted(model.getCompartment(c_id).getName() for c_id in get_r_comps(r.id, model)))),
                  get_kegg_r_id(r), get_gene_association(r), ','.join(get_pathway_expression(r)))
        if blocked_reactions:
            record += 'blocked' if r.id in blocked_reactions else '',
        data.append(record)
        index.append(r.id)

    columns = ["Id", "Name", "Lower Bound", "Upper Bound", "Formula", 'Compartments', "KEGG", "Gene association",
               "Subsystems"]
    if blocked_reactions:
        columns += ['Is blocked']
    return DataFrame(data=data, index=index, columns=columns)


def serialize_common_elements_to_csv(model_id2dfs, model_id2c_id_groups, model_id2m_id_groups,
                                     model_id2r_id_groups, prefix):
    def serialize_common_subpart_to_csv(i, model_id2id_groups, suffix):
        data = []

        for model_id2ids in model_id2id_groups:
            for model_id in sorted(model_id2ids.keys()):
                df = model_id2dfs[model_id][i]
                for el_id in sorted(model_id2ids[model_id]):
                    data_entry = [model_id]
                    data_entry.extend(df[df.Id == el_id].values[0])
                    data.append(data_entry)
            data.append([None])

        csv = '%s%s.csv' % (prefix, suffix)
        columns = ['Model']
        columns.extend(next(iter(model_id2dfs.values()))[i].columns)
        df2csv(DataFrame(data=data, columns=columns), csv)

        return csv

    m_csv = serialize_common_subpart_to_csv(0, model_id2m_id_groups, 'metabolites') \
        if model_id2m_id_groups else (0, None)
    comp_csv = serialize_common_subpart_to_csv(2, model_id2c_id_groups, 'compartments') \
        if model_id2c_id_groups else (0, None)
    r_csv = serialize_common_subpart_to_csv(1, model_id2r_id_groups, 'reactions') \
        if model_id2r_id_groups else (0, None)

    return comp_csv, m_csv, r_csv
