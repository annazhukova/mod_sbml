import getopt
import logging
import os
from os.path import dirname, abspath, basename
import sys

import libsbml

from em.em_manager import perform_efma

LINE_END = ' .'

ENTITY_IDENTIFIER_DELIMITER = ' '

REACTANT_PRODUCT_DELIMITER = ' = '

R_ID_FORMULA_DELIMITER = ' : '

METABOLITE_STOICHIOMETRY_DELIMITER = ' '

help_message = '''
Converts a model in dat format into SBML and (optionally) performs EFM analysis on it.
Specify the path to your model as the --dat parameter;
If you want to perform the EFM analysis, specify a reaction that should be present in EFMs as the --reaction parameter
and specify the path to your TreeEFM installation as the --tree parameter.
You can (optionally) specify the minimal motif length as the --motif parameter.
Usage: dat2sbml.py --dat model.dat --reaction AA2 --verbose --motif 3 --tree TreeEFMseq
'''


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

__author__ = 'anna'

R_REV = "-ENZREV"
R_IRREV = "-ENZIRREV"
M_INT = "-METINT"
M_EXT = "-METEXT"
R_DESCR = "-CAT"

HEADERS = [R_REV, R_IRREV, M_INT, M_EXT, R_DESCR]


def convert_metabolite(m_id, c_id, model, boundary=False, create_boundary_reaction=False):
    if m_id:
        m_name = m_id
        m_id = normalize_id(m_id)
        s = model.createSpecies()
        s.setId(m_id)
        s.setName(m_name)
        s.setCompartment(c_id)
        if boundary:
            if create_boundary_reaction:
                s = model.createSpecies()
                s.setId(m_id + '_b')
                s.setName(m_name)
                s.setCompartment(c_id)
                convert_reaction(model, 'r_%s_%s' % (m_id, s.getId()), True, ((m_id, 1),), ((s.getId(), 1),))
            s.setBoundaryCondition(boundary)


def normalize_id(id_, prefix='s_'):
    id_ = ''.join(e for e in id_ if e.isalnum())
    if id_[0].isdigit():
        id_ = prefix + id_
    return id_


def add_metabolites_to_reaction(r, m_id2st, is_reactant=True):
    for m_id, st in m_id2st:
        sr = r.createReactant() if is_reactant else r.createProduct()
        sr.setSpecies(m_id)
        sr.setStoichiometry(st)


def extract_m_id_stoichiometry_pairs(ms):
    for m_id in ms.split(' + '):
        m_id = m_id.strip()
        if m_id.find(METABOLITE_STOICHIOMETRY_DELIMITER) != -1:
            st, m_id = m_id.split(METABOLITE_STOICHIOMETRY_DELIMITER)
            st = float(st)
        else:
            st = 1
        m_id = normalize_id(m_id)
        yield m_id, st


def convert_reaction(model, r_id, rev, r_m_id2st, p_m_id2st):
    if r_id and (r_m_id2st or p_m_id2st):
        r = model.createReaction()
        r_name = r_id
        r_id = normalize_id(r_id, 'r_')
        r.setId(r_id)
        r.setName(r_name)
        r.setReversible(rev)
        add_metabolites_to_reaction(r, r_m_id2st, is_reactant=True)
        add_metabolites_to_reaction(r, p_m_id2st, is_reactant=False)


def convert_dat2sbml(in_dat, out_sbml, create_boundary_reaction=True):
    file_name, _ = os.path.splitext(os.path.basename(in_dat))
    file_name = normalize_id(file_name, 'M_') if file_name else 'Model'
    document = libsbml.SBMLDocument(2, 4)
    model = document.createModel()
    model.setId(file_name)
    c = model.createCompartment()
    c.setId('c')
    mode = None
    r_ids_rev = set()
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
                r_ids_rev = set(line.split(ENTITY_IDENTIFIER_DELIMITER))
            elif R_IRREV == mode:
                continue
            elif M_INT == mode:
                for m_id in set(line.split(ENTITY_IDENTIFIER_DELIMITER)):
                    convert_metabolite(m_id, c.getId(), model)
            elif M_EXT == mode:
                for m_id in set(line.split(ENTITY_IDENTIFIER_DELIMITER)):
                    convert_metabolite(m_id, c.getId(), model, True, create_boundary_reaction=create_boundary_reaction)
            elif R_DESCR == mode:
                line = line.replace(LINE_END, '')
                r_id, formula = line.split(R_ID_FORMULA_DELIMITER)
                rs, ps = formula.split(REACTANT_PRODUCT_DELIMITER)
                convert_reaction(model, r_id, r_id in r_ids_rev, extract_m_id_stoichiometry_pairs(rs),
                                 extract_m_id_stoichiometry_pairs(ps))
    libsbml.SBMLWriter().writeSBMLToFile(document, out_sbml)


def process_args(argv):
    try:
        opts, args = getopt.getopt(argv[1:], "m:r:h:v:d:t",
                                   ["help", "motif=", "reaction=", "verbose", "dat=", "tree="])
    except getopt.error, msg:
        raise Usage(msg)
    dat, r_id, verbose, motif_len, tree = None, None, False, 2, None
    # option processing
    for option, value in opts:
        if option in ("-h", "--help"):
            raise Usage(help_message)
        if option in ("-d", "--dat"):
            dat = value
        if option in ("-r", "--reaction"):
            r_id = value
        if option in ("-v", "--verbose"):
            verbose = True
        if option in ("-m", "--motif"):
            motif_len = int(value)
        if option in ("-t", "--tree"):
            tree = value
    if not dat or (r_id and not tree):
        raise Usage(help_message)
    model_dir = dirname(abspath(dat))
    sbml = '%s/%s.xml' % (model_dir, os.path.splitext(basename(dat))[0])
    return dat, r_id, sbml, model_dir, verbose, motif_len, tree


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        dat, r_id, sbml, model_dir, verbose, motif_len, tree = process_args(argv)
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2
    if verbose:
        logging.basicConfig(level=logging.INFO, format='%(message)s')
    convert_dat2sbml(dat, sbml)
    if r_id:
        perform_efma(sbml=sbml, in_r_id=r_id, in_r_reversed=False, out_r_id2rev_2threshold=(({r_id: False}, 1.0), ),
                     em_number=10000, min_motif_len=motif_len, model_dir=model_dir,
                     output_motif_file='%s/motifs.xlsx' % model_dir, output_efm_file='%s/efms.xlsx' % model_dir,
                     tree_efm_path=tree)


if __name__ == "__main__":
    sys.exit(main())