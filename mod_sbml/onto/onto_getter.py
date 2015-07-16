__author__ = 'anna'

import os

from mod_sbml import onto


def get_chebi():
    return os.path.join(os.path.dirname(os.path.abspath(onto.__file__)), '..', 'data', 'chebi.obo')


def get_go():
    return os.path.join(os.path.dirname(os.path.abspath(onto.__file__)), '..', 'data', 'gene_ontology_ext.obo')