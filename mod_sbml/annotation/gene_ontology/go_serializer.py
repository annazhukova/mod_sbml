import os

from mod_sbml.onto import parse, save_simple
from mod_sbml.onto.obo_ontology import PART_OF
from mod_sbml import annotation

__author__ = 'anna'

GO = os.path.join(os.path.dirname(os.path.abspath(annotation.__file__)), '..', 'data', 'go.obo')
SUB_GO_TXT = os.path.join(os.path.dirname(os.path.abspath(annotation.__file__)), '..', 'data', 'subgo.txt')
SUB_GO = os.path.join(os.path.dirname(os.path.abspath(annotation.__file__)), '..', 'data', 'subgo.obo')

CELLULAR_COMPONENT = 'go:0005575'


def get_go():
    return SUB_GO_TXT


def serialize_sub_go():
    onto = parse(GO, {PART_OF})
    terms_to_keep = onto.get_sub_tree(onto.get_term(CELLULAR_COMPONENT), relationships={PART_OF})
    for term in (t for t in onto.get_all_terms() if t not in terms_to_keep):
        onto.remove_term(term, True)
    save_simple(onto, SUB_GO_TXT)
    # save(onto, SUB_GO)


# serialize_sub_go()
