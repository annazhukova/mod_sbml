import getopt
import logging
import os
from os.path import dirname, abspath, basename
import sys
from constraint_based_analysis.efm.efm_analyser import perform_efma
from constraint_based_analysis.efm.efm_manager import efm2binary, get_int_size

from constraint_based_analysis.efm.metatool_manager import convert_metatool_output2efm, convert_dat2sbml

help_message = '''
Converts a model in dat format into SBML and (optionally) performs EFM analysis on it.
Specify the path to your model as the --dat parameter;

If you want to perform the EFM calculation followed by the EFM classification,
specify a reaction that should be present in EFMs as the --reaction parameter
and specify the path to your TreeEFM installation as the --tree parameter.

If you want to perform a classification of your EFMs obtained with Metatool,
specify the path to your Metatool output file as the --efm parameter.
If you are only interested in EFMs containing a particular reaction, specify its id as the --reaction parameter.

Usage:  (Conversion to SBML + Calculation of EFMs + EFM classification)
            main.py --dat model.dat --reaction AA2 --verbose --tree TreeEFMseq
        (Conversion to SBML + Classification of EFMs obtained with Metatool)
            main.py --dat model.dat --reaction AA2 --verbose --efm out-meta
'''


class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

__author__ = 'anna'


def process_args(argv):
    try:
        opts, args = getopt.getopt(argv[1:], "r:h:v:d:t:e",
                                   ["help", "reaction=", "verbose", "dat=", "tree=", "efm="])
    except getopt.error, msg:
        raise Usage(msg)
    dat, r_id, verbose, tree, efm = None, None, False, None, None
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
        if option in ("-t", "--tree"):
            tree = value
        if option in ("-e", "--efm"):
            efm = value
    if not dat or (not efm and (not tree or (tree and not r_id))):
        raise Usage(help_message)
    model_dir = dirname(abspath(dat))
    sbml = '%s/%s.xml' % (model_dir, os.path.splitext(basename(dat))[0])
    return dat, r_id, sbml, model_dir, verbose, tree, efm

def main(argv=None):
    if argv is None:
        argv = sys.argv
    try:
        dat, r_id, sbml, model_dir, verbose, tree, efm_file = process_args(argv)
    except Usage, err:
        print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
        print >> sys.stderr, "\t for help use --help"
        return 2
    if verbose:
        logging.basicConfig(level=logging.INFO, format='%(asctime)s: %(message)s')

    r_rev_ids, r_irrev_ids = convert_dat2sbml(dat, sbml, create_boundary_reaction=True)
    r_ids = r_rev_ids + r_irrev_ids
    if efm_file:
        efms = convert_metatool_output2efm(efm_file, r_rev_ids=r_rev_ids, r_irrev_ids=r_irrev_ids)
        if r_id:
            efms = [efm for efm in efms if r_id in efm]
        int_size = get_int_size()
        efms = [efm2binary(efm, r_ids, r_rev_ids, int_size) for efm in efms]
    else:
        efms = None

    perform_efma(in_r_id=r_id, in_r_reversed=False, out_r_id2rev_2threshold=(), sbml=sbml,
                 directory=model_dir, tree_efm_path=tree, efms=efms, r_ids=r_ids,
                 output_efm_file='%s/efms.txt' % model_dir, imp_rn_threshold=0,
                 calculate_important_reactions=True,
                 calculate_patterns=True)


if __name__ == "__main__":
    sys.exit(main())