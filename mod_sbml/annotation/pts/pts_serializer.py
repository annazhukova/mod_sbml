import os

from pandas import notnull

from mod_sbml import annotation
from mod_sbml.serialization import csv2df
from mod_sbml.onto import Ontology, Term, save_simple
from mod_sbml.onto.obo_ontology import PART_OF

DELIMITER = '|'

__author__ = 'anna'

PTS_CSV = os.path.join(os.path.dirname(os.path.abspath(annotation.__file__)), '..', 'data', 'PTS.csv')
SUB_PTS_TXT = os.path.join(os.path.dirname(os.path.abspath(annotation.__file__)), '..', 'data', 'PTS.txt')

METABOLIC_PATHWAY = 'http://purl.obolibrary.org/obo/iev_0000818'


def get_pts():
    return SUB_PTS_TXT


def parse_csv(csv_file=PTS_CSV, sep=','):
    if not csv_file or not os.path.exists(csv_file):
        raise ValueError('csv file was not found in %s' % csv_file)
    df = csv2df(csv_file, sep=sep)
    df = df.where((notnull(df)), None)
    ontology = Ontology()

    for t_id, row in df.iterrows():
        name, synonyms, parents, wholes, xrefs = row['Preferred Label'], row['http://scai.fraunhofer.de/PTS#Synonym'], \
                                                 row['Parents'], row['http://data.bioontology.org/metadata/obo/part_of'], \
                                                 row['http://purl.obolibrary.org/obo/xref_analog']
        if not t_id:
            continue
        term = Term(onto=ontology)
        term.set_id(t_id)
        if name:
            term.set_name(name.replace('\\"', '"'))
        if parents:
            for parent in parents.split(DELIMITER):
                term.add_parent(parent)
        if wholes:
            for whole in wholes.split(DELIMITER):
                ontology.add_relationship(term.get_id(), PART_OF, whole)
        if synonyms:
            for synonym in synonyms.split(DELIMITER):
                term.add_synonym(synonym)
        if xrefs:
            for xref in xrefs.split(DELIMITER):
                last_colon_index = xref.rfind(':')
                if last_colon_index == -1 or last_colon_index == len(xref):
                    continue
                db_name, value = xref[:last_colon_index], xref[last_colon_index + 1:]
                term.add_xref(db_name, value)
        ontology.add_term(term)
    return ontology


def serialize_sub_pts():
    onto = parse_csv()
    terms_to_keep = onto.get_descendants(METABOLIC_PATHWAY, False)
    for term in (t for t in onto.get_all_terms() if t.get_id() not in terms_to_keep):
        onto.remove_term(term, True)
    save_simple(onto, SUB_PTS_TXT)


# serialize_sub_pts()
