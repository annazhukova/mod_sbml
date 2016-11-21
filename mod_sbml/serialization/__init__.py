from pandas import DataFrame
from mod_sbml.sbml.sbml_manager import get_reactants, get_products

__author__ = 'anna'


def df2csv(df, path, sep='\t'):
    df.to_csv(path_or_buf=path, na_rep='', sep=sep, index=False)


def csv2df(path, sep='\t'):
    return DataFrame.from_csv(path, sep=sep, index_col=0)


def get_sbml_r_formula(model, r, show_compartments=True, show_metabolite_ids=True):
    format_m = lambda m_id, st: "%s%s" % (
        "%.2g " % st if st != 1 else "",
        format_m_name(model.getSpecies(m_id), model, show_compartments, show_metabolite_ids))
    formula = " + ".join([format_m(m_id, st) for (m_id, st) in get_reactants(r, True)]) + \
              (" <=> " if r.getReversible() else " --> ") + \
              " + ".join([format_m(m_id, st) for (m_id, st) in get_products(r, True)])
    return formula


def format_m_name(m, model, show_compartment=True, show_id=True):
    name = m.getName()
    c = model.getCompartment(m.getCompartment())
    c_name = c.getName()
    if -1 == name.find(c_name) and show_compartment:
        name = "%s[%s]" % (name, c_name)
    if not show_compartment:
        name = name.replace("[%s]" % c_name, "").replace("[%s]" % c.getId(), "").strip()
    return "%s(%s)" % (name, m.id) if show_id else name