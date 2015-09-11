from collections import defaultdict

import libsbml
from pandas import DataFrame

from kegg.kegg_annotator import get_kegg_r_id, get_kegg_m_id
from reaction_boundary_manager import get_bounds
from sbml_manager import get_gene_association
from serialization_manager import get_sbml_r_formula

__author__ = 'anna'


def serialize_model_info(sbml, prefix):
    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()

    def to_csv(to_df, name):
        csv = '%s%s.csv' % (prefix, name)
        df2csv(to_df(model), csv)
        return csv

    return to_csv(compartments2df, 'compartments'), to_csv(metabolites2df, 'metabolites'), \
           to_csv(reactions2df, 'reactions')


def df2csv(df, path):
    df.to_csv(path_or_buf=path, na_rep='', sep='\t', index=False)


def metabolites2df(model):
    data = []
    index = []
    for m in sorted(model.getListOfSpecies(), key=lambda m: m.id):
        c = model.getCompartment(m.getCompartment())
        data.append((m.id, m.name, c.name if c.name else c.id, get_kegg_m_id(m)))
        index.append(m.id)
    return DataFrame(data=data, index=index, columns=['Id', 'Name', 'Compartment', 'KEGG'])


def compartments2df(model):
    data = []
    index = []
    for c in sorted(model.getListOfCompartments(), key=lambda c: c.id):
        data.append((c.id, c.name if c.name else c.id,))
        index.append(c.id)
    return DataFrame(data=data, index=index, columns=['Id', 'Name'])


def reactions2df(model):
    data = []
    index = []
    for r in sorted(model.getListOfReactions(), key=lambda r: r.id):
        data.append((r.id, r.name, get_bounds(r)[0], get_bounds(r)[1], get_sbml_r_formula(model, r, False),
                     get_kegg_r_id(r), get_gene_association(r)))
        index.append(r.id)
    return DataFrame(data=data, index=index,
                     columns=["Id", "Name", "Lower Bound", "Upper Bound", "Formula", "Kegg", "Gene association"])


def serialize_common_part_to_csv(merged_sbml, sbml2id2id, common_ids, sbml2name, prefix):
    doc = libsbml.SBMLReader().readSBML(merged_sbml)
    merged_model = doc.getModel()
    id_sbml2id = defaultdict(list)
    for sbml, id2id in sbml2id2id.iteritems():
        for (s_id, t_id) in id2id.iteritems():
            id_sbml2id[(t_id, sbml)].append(s_id)

    def serialize_common_subpart_to_csv(get_df, suffix):
        data = []
        num = 0
        sbml2df = {}

        for sbml in sbml2id2id.iterkeys():
            d = libsbml.SBMLReader().readSBML(sbml)
            model = d.getModel()
            sbml2df[sbml] = get_df(model)

        df = get_df(merged_model)

        for element in df.values:
            element_id = element[0]
            if element_id in common_ids:
                num += 1
                data_entry = ['Merged model']
                data_entry.extend(element)
                data.append(data_entry)
                for sbml in sbml2id2id.iterkeys():
                    if (element_id, sbml) in id_sbml2id:
                        df = sbml2df[sbml]
                        for id_ in id_sbml2id[(element_id, sbml)]:
                            data_entry = [sbml2name[sbml]]
                            data_entry.extend(df[df.Id == id_].values[0])
                            data.append(data_entry)
                data.append([None])

        csv = '%s%s.csv' % (prefix, suffix)
        columns = ['Model']
        columns.extend(df.columns)
        df2csv(DataFrame(data=data, columns=columns), csv)

        return num, csv

    cc_num, comp_csv = serialize_common_subpart_to_csv(compartments2df, 'compartments')
    cm_num, m_csv = serialize_common_subpart_to_csv(metabolites2df, 'metabolites')
    cr_num, r_csv = serialize_common_subpart_to_csv(reactions2df, 'reactions')

    return (cc_num, cm_num, cr_num), (comp_csv, m_csv, r_csv)
