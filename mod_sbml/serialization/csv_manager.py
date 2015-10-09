from pandas import DataFrame

from mod_sbml.annotation.kegg.kegg_annotator import get_kegg_r_id, get_kegg_m_id
from mod_sbml.sbml.reaction_boundary_manager import get_bounds
from mod_sbml.sbml.sbml_manager import get_gene_association, get_formulas
from mod_sbml.serialization import get_sbml_r_formula, df2csv

__author__ = 'anna'


def serialize_model_info(model, prefix):

    def to_csv(to_df, name):
        csv = '%s%s.csv' % (prefix, name)
        df2csv(to_df(model), csv)
        return csv

    return to_csv(compartments2df, 'compartments'), to_csv(metabolites2df, 'metabolites'), \
           to_csv(reactions2df, 'reactions')


def metabolites2df(model, m_id2chebi_id=None):
    data = []
    index = []
    for m in sorted(model.getListOfSpecies(), key=lambda m: m.id):
        c = model.getCompartment(m.getCompartment())
        formulas = get_formulas(m)
        data_entry = (m.id, m.name, c.id, get_kegg_m_id(m), formulas.pop() if formulas else None)
        if m_id2chebi_id:
            data_entry += (m_id2chebi_id[m.id] if m.id in m_id2chebi_id else None,)
        data.append(data_entry)
        index.append(m.id)
    columns = ['Id', 'Name', 'Compartment', 'KEGG', 'Formula']
    if m_id2chebi_id:
        columns.append('ChEBI')
    return DataFrame(data=data, index=index, columns=columns)


def compartments2df(model, get_term=None):
    data = []
    index = []
    for c in sorted(model.getListOfCompartments(), key=lambda c: c.id):
        data_piece = (c.id, c.name, get_term(c)) if get_term else (c.id, c.name)
        data.append(data_piece)
        index.append(c.id)
    return DataFrame(data=data, index=index, columns=['Id', 'Name', "GO"] if get_term else ['Id', 'Name'])


def reactions2df(model):
    data = []
    index = []
    for r in sorted(model.getListOfReactions(), key=lambda r: r.id):
        data.append((r.id, r.name, get_bounds(r)[0], get_bounds(r)[1], get_sbml_r_formula(model, r, False),
                     get_kegg_r_id(r), get_gene_association(r)))
        index.append(r.id)
    return DataFrame(data=data, index=index,
                     columns=["Id", "Name", "Lower Bound", "Upper Bound", "Formula", "Kegg", "Gene association"])


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
        columns.extend(next(model_id2dfs.itervalues())[i].columns)
        df2csv(DataFrame(data=data, columns=columns), csv)

        return csv

    m_csv = serialize_common_subpart_to_csv(0, model_id2m_id_groups, 'metabolites') \
        if model_id2m_id_groups else (0, None)
    comp_csv = serialize_common_subpart_to_csv(2, model_id2c_id_groups, 'compartments') \
        if model_id2c_id_groups else (0, None)
    r_csv = serialize_common_subpart_to_csv(1, model_id2r_id_groups, 'reactions') \
        if model_id2r_id_groups else (0, None)

    return comp_csv, m_csv, r_csv
