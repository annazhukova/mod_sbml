import os

from mod_sbml.onto.obo_ontology import Ontology
from mod_sbml.onto.term import Term, FORMULA

__author__ = 'anna'

RELS_HEADER = '[Relationships]\n'
TERMS_HEADER = '[Terms]\n'


def parse_simple(path):
    if not os.path.exists(path):
        return None
    ontology = Ontology()
    TERMS = 0
    RELS = 1
    mode = None
    with open(path, 'r') as f:
        for line in f:
            if not line or '\n' == line.strip():
                continue
            if line == TERMS_HEADER:
                mode = TERMS
                continue
            if line == RELS_HEADER:
                mode = RELS
                continue
            if TERMS == mode:
                term = Term(onto=ontology, s=line)
                ontology.add_term(term)
            elif RELS == mode:
                line = line.splitlines()[0]
                subj, rel, obj = line.split('\t')
                ontology.add_relationship(subj=subj, rel=rel, obj=obj)
    return ontology


def parse(obo_file, relationships=None):
    if not obo_file or obo_file.find(".obo") == -1 or not os.path.exists(obo_file):
        return None
    ontology = Ontology()
    term = None
    with open(obo_file, 'r') as obo:
        for line in obo:
            line = line.replace("\n", '')
            if line.find("[Term]") != -1:
                if term:
                    ontology.add_term(term)
                term = Term(onto=ontology)
                continue
            if not term:
                continue
            semicolon = line.find(":")
            if semicolon == -1:
                continue
            prefix = line[0:semicolon].strip()
            value = line[semicolon + 1: len(line)].strip()
            comment = value.find("!")
            if comment != -1:
                value = value[0:comment].strip()
            if prefix == "id":
                term.set_id(value)
            elif prefix == "is_obsolete" and value == "true":
                term = None
            elif prefix == "alt_id":
                term.add_alt_id(value)
            elif prefix == "name":
                term.set_name(value.replace('\\"', '"'))
            elif prefix == "is_a":
                term.add_parent(value)
            elif prefix == "relationship":
                comment_start = value.find("!")
                if comment_start != -1:
                    value = value[0:comment_start]
                first, second = value.strip().split(" ")
                if not relationships or first in relationships:
                    ontology.add_relationship(term.get_id(), first, second.lower())
            elif prefix == "synonym":
                start = value.find('"')
                if start == -1:
                    continue
                st = start + 1
                while -1 != value.find('\\"', st):
                    st = value.find('\\"', st) + 2
                end = value.find('"', st)
                if end == -1:
                    continue
                synonym = value[start + 1:end].replace('\\"', '"')
                if synonym and synonym != '.':
                    term.add_synonym(synonym)
                    if end + 1 < len(value) and -1 != value[end + 1:].find('FORMULA'):
                        term.add_xref(FORMULA, synonym)
            elif prefix == "xref":
                value = value.strip()
                if not value:
                    continue
                comment = value.find('"')
                if comment != - 1:
                    value = value[:comment]
                modifier = value.find("{")
                if modifier != -1:
                    value = value[:modifier]
                colon = value.find(':')
                if -1 != colon:
                    db_name = value[:colon]
                    value = value[colon + 1:].strip()
                    if db_name and value:
                        term.add_xref(db_name, value)
    if term:
        ontology.add_term(term)
    return ontology


def filter_ontology(onto, terms_collection, relationships=None, min_deepness=None):
    terms_to_keep = set()

    for term in terms_collection:
        if term in terms_to_keep:
            continue
        terms_to_keep |= onto.get_sub_tree(term, relationships=relationships)
        for ancestor in onto.get_generalized_ancestors(term, direct=False, checked=set(), relationships=relationships,
                                                       depth=min_deepness):
            terms_to_keep |= onto.get_sub_tree(ancestor, relationships=relationships)

    for term in (t for t in onto.get_all_terms() if t not in terms_to_keep):
        onto.remove_term(term, True)

    onto.filter_relationships(relationships)


def save_simple(onto, path):
    processed_terms = set()
    with open(path, 'w') as f:
        f.write(TERMS_HEADER)
        for term in onto.get_all_terms():
            id_ = term.get_id()
            if id_ in processed_terms:
                continue
            processed_terms.add(id_)
            f.write(str(term))
        f.write(RELS_HEADER)
        for rel in reduce(lambda s1, s2: s1 | s2, onto.rel_map.itervalues(), set()):
            f.write("%s\t%s\t%s\n" % rel)


def save(onto, path):
    processed = set()
    with open(path, 'w') as f:
        f.write('format-version: 1.2\n')
        f.write('synonymtypedef: UNKNOWN "Unknown"\n')
        f.write('\n')
        for term in onto.get_all_terms():
            id_ = term.get_id()
            if id_ in processed:
                continue
            processed.add(id_)
            f.write("[Term]\n")
            f.write("id: {0}\n".format(id_))
            for a_id in term.get_all_ids():
                if a_id == id_:
                    continue
                f.write("alt_id: {0}\n".format(a_id))
            name = term.get_name()
            f.write("name: {0}\n".format(name))
            for syn in term.get_synonyms():
                if syn == name:
                    continue
                f.write('synonym: "{0}" RELATED UNKNOWN [ChEBI:]\n'.format(syn.replace('"', '\\"')))
            for parent in term.get_parent_ids():
                f.write("is_a: {0}\n".format(parent))
            for db in term.get_dbs():
                for value in term.get_xrefs(db):
                    f.write('xref: %s:%s "%s"' % (db, value, db))
            for (subj, rel, obj) in onto.get_term_relationships(id_, None, 1):
                f.write("relationship: {1} {0}\n".format(obj, rel))
            f.write("\n")
        for rel in onto.get_relationships():
            f.write("[Typedef]\n")
            f.write("id: {0}\n".format(rel))
            f.write("name: {0}\n".format(rel.replace("_", " ")))
            f.write("is_cyclic: false\n")
            f.write("is_transitive: false\n")