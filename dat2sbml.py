import getopt
import logging
import os
from os.path import dirname, abspath, basename
import sys

import libsbml
import re

from em.em_manager import perform_efma, classify_efms, serialize_efms

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
                convert_reaction(model, 'r_%s_in_out' % m_id, True, ((m_id, 1),), ())
            else:
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


def convert_metatool_output2efm(metatool_file, sorted_r_ids):
    """
    Extracts the EFM information from a given output file of Metatool
    and converts them to an inner format: list of EFMs represented as dictionaries of r_id: coefficient,
    i.e. [{r_id: coefficient, ...}, ...]
    :param metatool_file: output file produced by Metatool
    :param sorted_r_ids: list of reaction ids, in the same order as in the input .dat file given to Metatool,
    reversible reactions followed by irreversible.
    :return: list of EFMs represented as dictionaries of r_id: coefficient,
    i.e. [{r_id: coefficient, ...}, ...]
    """
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
                efm = {r_id: float(c) for (r_id, c) in zip(sorted_r_ids, coefficients) if float(c) != 0}
                efms.append(efm)
            line = next(f)
    return efms


def convert_dat2sbml(in_dat, out_sbml, create_boundary_reaction=True):
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
    file_name = normalize_id(file_name, 'M_') if file_name else 'Model'
    document = libsbml.SBMLDocument(2, 4)
    model = document.createModel()
    model.setId(file_name)
    c = model.createCompartment()
    c.setId('c')
    mode = None
    r_ids_rev = set()
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
                r_rev_ids = line.split(ENTITY_IDENTIFIER_DELIMITER)
                r_ids_rev = set(r_rev_ids)
            elif R_IRREV == mode:
                r_irrev_ids = line.split(ENTITY_IDENTIFIER_DELIMITER)
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
    logging.info("...%s converted to %s" % (in_dat, out_sbml))
    return r_rev_ids + r_irrev_ids


def process_args(argv):
    try:
        opts, args = getopt.getopt(argv[1:], "m:r:h:v:d:t:e",
                                   ["help", "motif=", "reaction=", "verbose", "dat=", "tree=", "efm="])
    except getopt.error, msg:
        raise Usage(msg)
    dat, r_id, verbose, motif_len, tree, efm = None, None, False, 2, None, None
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
        if option in ("-e", "--efm"):
            efm = value
    if not dat:
        raise Usage(help_message)
    model_dir = dirname(abspath(dat))
    sbml = '%s/%s.xml' % (model_dir, os.path.splitext(basename(dat))[0])
    return dat, r_id, sbml, model_dir, verbose, motif_len, tree, efm


def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        dat, r_id, sbml, model_dir, verbose, motif_len, tree, efm_file = process_args(argv)
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2
    if verbose:
        logging.basicConfig(level=logging.INFO, format='%(message)s')
    r_ids = convert_dat2sbml(dat, sbml, create_boundary_reaction=False)
    if tree:
        perform_efma(sbml=sbml, in_r_id=r_id, in_r_reversed=False, out_r_id2rev_2threshold=(({r_id: False}, 1.0), ),
                     em_number=10000, min_motif_len=motif_len, model_dir=model_dir,
                     output_motif_file='%s/motifs.xlsx' % model_dir, output_efm_file='%s/efms.xlsx' % model_dir,
                     tree_efm_path=tree)
    if efm_file:
        efms = convert_metatool_output2efm(efm_file, r_ids)
        if r_id:
            efms = [efm for efm in efms if r_id in efm]
        serialize_efms(sbml, efms, path='%s/efms.xlsx' % model_dir)
        clusters, outliers = classify_efms(efms, min_motif_length=motif_len,
                                           r_ids=r_ids, output_file='%s/motifs.xlsx' % model_dir, sbml=sbml)


if __name__ == "__main__":
    sys.exit(main())