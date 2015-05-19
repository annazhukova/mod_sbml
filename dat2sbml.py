import getopt
import logging
import os
from os.path import dirname, abspath, basename
import sys

import libsbml

from em.em_manager import perform_efma


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


def convert_metabolite(m_id, c_id, model, boundary=False):
    if m_id:
        if m_id[0].isdigit():
            m_id = 's_' + m_id
        s = model.createSpecies()
        s.setId(m_id)
        s.setCompartment(c_id)
        s.setBoundaryCondition(boundary)


def normalize_id(id_, prefix='s_'):
    if id_[0].isdigit():
        id_ = prefix + id_
    return id_


def add_metabolites_to_reaction(r, ms, is_reactant=True):
    for m_id in ms:
        m_id = m_id.strip()
        if m_id.find(' ') != -1:
            st, m_id = m_id.split(' ')
            st = float(st)
        else:
            st = 1
        m_id = normalize_id(m_id)
        sr = r.createReactant() if is_reactant else r.createProduct()
        sr.setSpecies(m_id)
        sr.setStoichiometry(st)


def convert_reaction(model, r_id, rev, rs, ps):
    r = model.createReaction()
    r_id = normalize_id(r_id, 'r_')
    r.setId(r_id)
    r.setReversible(rev)
    add_metabolites_to_reaction(r, rs, is_reactant=True)
    add_metabolites_to_reaction(r, ps, is_reactant=False)


def convert_dat2sbml(in_dat, out_sbml):
    file_name, file_extension = os.path.splitext(os.path.basename(in_dat))
    if file_name and file_name[0].isdigit():
        file_name = 'M_' + file_name
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
                r_ids_rev = set(line.split(' '))
            elif R_IRREV == mode:
                continue
            elif M_INT == mode:
                for m_id in set(line.split(' ')):
                    convert_metabolite(m_id, c.getId(), model)
            elif M_EXT == mode:
                for m_id in set(line.split(' ')):
                    convert_metabolite(m_id, c.getId(), model, True)
            elif R_DESCR == mode:
                line = line.replace(' .', '')
                r_id, formula = line.split(' : ')
                rs, ps = formula.split(' = ')
                rs, ps = rs.split(' + '), ps.split(' + ')
                convert_reaction(model, r_id, r_id in r_ids_rev, rs, ps)
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