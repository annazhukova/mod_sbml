from collections import Counter, defaultdict
import logging
import os
import re

import libsbml

from mod_sbml.constraint_based_analysis.efm.EFM import EFM
from mod_sbml.gibbs.reaction_boundary_manager import set_bounds
from serialization_manager import get_sbml_r_formula

__author__ = 'anna'

LINE_END = ' .'
ENTITY_IDENTIFIER_DELIMITER = ' '
REACTANT_PRODUCT_DELIMITER = ' = '
R_ID_FORMULA_DELIMITER = ' : '
METABOLITE_STOICHIOMETRY_DELIMITER = ' '

R_REV = "-ENZREV"
R_IRREV = "-ENZIRREV"
M_INT = "-METINT"
M_EXT = "-METEXT"
R_DESCR = "-CAT"

HEADERS = [R_REV, R_IRREV, M_INT, M_EXT, R_DESCR]

def get_c_id(model, m_id):
    if model.getCompartment(m_id[-1:]):
        return m_id[-1:]
    return None

def convert_metabolite(m_id, model, boundary=False, create_boundary_reaction=False, abbr2name=None, c_id=None):
    if m_id:
        if not c_id and model.getCompartment(m_id[-1:]):
            c_id = [m_id[-1:]]
            m_name = m_id[:-1]
        else:
            m_name = m_id
        if abbr2name and m_name in abbr2name:
            m_name = abbr2name[m_name]
        m_id = normalize_id(m_id)
        s = model.createSpecies()
        s.setId(m_id)
        s.setName(m_name)
        s.setCompartment(c_id)
        if boundary:
            if create_boundary_reaction:
                convert_reaction(model, 'r_%s_exchange' % m_id, True, {m_id: 1}, {},
                                 r_name='Exchange of %s' % m_name)
            else:
                s.setBoundaryCondition(boundary)


def convert_metabolites(m_id2c_id_b, model, create_boundary_reaction=False, abbr2name=None):
    for m_id, (c_id, b) in m_id2c_id_b.iteritems():
        if m_id and c_id:
            m_name = m_id[:-1]
            if abbr2name and m_name in abbr2name:
                m_name = abbr2name[m_name]
            m_id_norm = normalize_id(m_id)
            s = model.createSpecies()
            s.setId(m_id_norm)
            s.setName(m_name)
            s.setCompartment(c_id)
            if b:
                if create_boundary_reaction:
                    convert_reaction(model, 'r_%s_exchange' % m_id_norm, True, {m_id: 1}, {},
                                     r_name='Exchange of %s' % m_name, m_id2c_id_b=m_id2c_id_b)
                else:
                    s.setBoundaryCondition(b)


def normalize_id(id_, prefix='s_'):
    id_ = ''.join(e for e in id_ if e.isalnum() or '_' == e)
    if id_ and id_[0].isdigit():
        id_ = prefix + id_
    return id_


def add_metabolites_to_reaction(r, m_id2st, is_reactant=True, default_c_id=None, m_id2c_id_b=None):
    for m_id, st in m_id2st.iteritems():
        sr = r.createReactant() if is_reactant else r.createProduct()
        if m_id2c_id_b:
            c_id, b = m_id2c_id_b[m_id]
            if c_id:
                sr.setSpecies(normalize_id(m_id))
            else:
                c_m_id = m_id + default_c_id
                sr.setSpecies(normalize_id(c_m_id))
                m_id2c_id_b[c_m_id] = default_c_id, b
        else:
            sr.setSpecies(normalize_id(m_id))
        sr.setStoichiometry(st)


def extract_m_id_stoichiometry_pairs(ms):
    for m_id in ms.split(' + '):
        m_id = m_id.strip()
        if m_id:
            if m_id.find(METABOLITE_STOICHIOMETRY_DELIMITER) != -1:
                st, m_id = m_id.split(METABOLITE_STOICHIOMETRY_DELIMITER)
                st = float(st)
            else:
                st = 1
            yield m_id, st


def convert_reaction(model, r_id, rev, r_m_id2st, p_m_id2st, default_c_id=None, r_name=None, m_id2c_id_b=None):
    if r_id and (r_m_id2st or p_m_id2st):
        r = model.createReaction()
        if not r_name:
            r_name = r_id
        r_id = normalize_id(r_id, 'r_')
        r.setId(r_id)
        r.setName(r_name)
        r.setReversible(rev)

        if m_id2c_id_b:
            comps = Counter([m_id2c_id_b[m_id][0] for m_id in r_m_id2st.iterkeys()
                             if m_id2c_id_b[m_id][0]] +
                            [m_id2c_id_b[m_id][0] for m_id in p_m_id2st.iterkeys()
                             if m_id2c_id_b[m_id][0]])
            if comps:
                default_c_id, _ = comps.most_common(1)[0]
        add_metabolites_to_reaction(r, r_m_id2st, is_reactant=True, default_c_id=default_c_id,
                                    m_id2c_id_b=m_id2c_id_b)
        add_metabolites_to_reaction(r, p_m_id2st, is_reactant=False, default_c_id=default_c_id,
                                    m_id2c_id_b=m_id2c_id_b)
        set_bounds(r, -10 if rev else 0, 10)


def convert_dat2sbml(in_dat, out_sbml, create_boundary_reaction=True, c_id2c_name=None, default_c_id=None,
                     abbr2name=None):
    """
    Converts a .dat file with the model (in the format accepted by Metatool) into SBML.
    :param in_dat: input .dat file
    :param out_sbml: path where to save the output SBML file
    :param create_boundary_reaction: set to False if you want the external metabolites to be marked as boundary;
    set to True if you want the external metabolites to considered as internal
    with fake input/output reactions (M ->; -> M) added;
    :return: list of reaction ids, in the same order as in the input .dat file,
    reversible reactions followed by irreversible.
    """
    logging.info("Converting %s to SBML..." % in_dat)
    file_name, _ = os.path.splitext(os.path.basename(in_dat))
    m_name = file_name if file_name else 'Model'
    m_id = normalize_id(m_name, 'M_')
    document = libsbml.SBMLDocument(2, 4)
    model = document.createModel()
    model.setId(m_id)
    model.setName(m_name)
    if not c_id2c_name:
        if default_c_id:
            c_id2c_name = {default_c_id: default_c_id}
        else:
            c_id2c_name = {'c': "cell"}
            default_c_id = 'c'
    for c_id, c_name in c_id2c_name.iteritems():
        c = model.createCompartment()
        c.setId(c_id)
        c.setName(c_name)
    mode = None
    r_ids_rev = set()
    r_rev_ids, r_irrev_ids = [], []
    m_id2c_id_b = {}
    with open(in_dat, 'r') as f:
        for line in f:
            line = line.replace('\n', '').strip()
            if not line:
                continue
            if line in HEADERS:
                mode = line
                continue
            if not mode:
                continue
            if R_REV == mode:
                r_rev_ids = [r_id for r_id in line.split(ENTITY_IDENTIFIER_DELIMITER) if r_id.strip()]
                r_ids_rev = set(r_rev_ids)
            elif R_IRREV == mode:
                r_irrev_ids = [r_id for r_id in line.split(ENTITY_IDENTIFIER_DELIMITER) if r_id.strip()]
            elif M_INT == mode:
                for m_id in set(line.split(ENTITY_IDENTIFIER_DELIMITER)):
                    if m_id:
                        m_id2c_id_b[m_id] = get_c_id(model, m_id), False
                    # convert_metabolite(m_id, model, abbr2name=abbr2name, m_id2c_ids=m_id2c_ids)
            elif M_EXT == mode:
                for m_id in set(line.split(ENTITY_IDENTIFIER_DELIMITER)):
                    if m_id:
                        m_id2c_id_b[m_id] = get_c_id(model, m_id), True
                    # convert_metabolite(m_id, model, True,
                    #                    create_boundary_reaction=create_boundary_reaction, abbr2name=abbr2name,
                    #                    m_id2c_ids=m_id2c_id_b)
            elif R_DESCR == mode:
                line = line.replace(LINE_END, '')
                r_id, formula = line.split(R_ID_FORMULA_DELIMITER)
                rs, ps = formula.split(REACTANT_PRODUCT_DELIMITER)
                convert_reaction(model, r_id, r_id in r_ids_rev, dict(extract_m_id_stoichiometry_pairs(rs)),
                                 dict(extract_m_id_stoichiometry_pairs(ps)), default_c_id=default_c_id,
                                 m_id2c_id_b=m_id2c_id_b)
    convert_metabolites(m_id2c_id_b, model, create_boundary_reaction=create_boundary_reaction, abbr2name=abbr2name)
    libsbml.SBMLWriter().writeSBMLToFile(document, out_sbml)
    logging.info("...%s converted to %s" % (in_dat, out_sbml))
    return r_rev_ids, r_irrev_ids


def convert_metatool_output2efm(metatool_file, in_dat=None, r_rev_ids=None, r_irrev_ids=None, r_id=None):
    """
    Extracts the EFM information from a given output file of Metatool
    and converts them to a list of EFMs.

    :param in_dat: input .dat file with the model (in the format accepted by Metatool),
    :param metatool_file: output file produced by Metatool,
    :param r_id: (optional) reaction id, if specified only the EFMs containing this reaction will be considered.
    :return: list of EFMs
    """
    if in_dat:
        r_rev_ids, r_irrev_ids = extract_r_ids_from_dat(in_dat)
    if not r_irrev_ids and not r_rev_ids:
        raise ValueError('You should either specify the reaction lists or the in_dat file!')
    sorted_r_ids = r_rev_ids + r_irrev_ids
    r_rev_ids = set(r_rev_ids)
    efms = []
    with open(metatool_file, 'r') as f:
        line = next(f)
        while -1 == line.find('ELEMENTARY MODES'):
            line = next(f)
        line = next(f)
        while -1 == line.find('matrix dimension'):
            line = next(f)
        line = next(f)
        while -1 == line.find('reactions (columns) are sorted in the same order as in the ENZREV ENZIRREV section'):
            line = line.replace('\n', '').strip()
            if line:
                coefficients = re.findall(r"[-+]?\d*\.*\d+", line)
                r_id2coefficient = {r_id: float(c) for (r_id, c) in zip(sorted_r_ids, coefficients) if float(c) != 0}
                if not r_id or r_id in r_id2coefficient:
                    efms.append(EFM(r_id2coeff=r_id2coefficient, r_ids=sorted_r_ids, rev_r_ids=r_rev_ids))
            line = next(f)
    return efms


def extract_r_ids_from_dat(in_dat):
    """
    Extracts a list of reversible and a list of irreversible reaction identifiers from a .dat file with the model
    (in the format accepted by Metatool).
    :param in_dat: input .dat file
    :return: list of reversible reaction ids, in the same order as in the input .dat file,
            list of irreversible reaction ids, in the same order as in the input .dat file.
    """
    file_name, _ = os.path.splitext(os.path.basename(in_dat))
    r_rev_ids, r_irrev_ids = [], []
    with open(in_dat, 'r') as f:
        for line in f:
            line = line.replace('\n', '').strip()
            if not line:
                continue
            if line in HEADERS:
                mode = line
                continue
            if not mode:
                continue
            if R_REV == mode:
                r_rev_ids = [r_id for r_id in line.split(ENTITY_IDENTIFIER_DELIMITER) if r_id.strip()]
            elif R_IRREV == mode:
                r_irrev_ids = [r_id for r_id in line.split(ENTITY_IDENTIFIER_DELIMITER) if r_id.strip()]
    return r_rev_ids, r_irrev_ids