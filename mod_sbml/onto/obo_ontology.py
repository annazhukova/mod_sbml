#!/usr/bin/env python
# encoding: utf-8
from collections import defaultdict

from mod_sbml.utils.misc import remove_from_map

__author__ = 'anna'

import os


def filter_ontology(onto, terms_collection, relationships=None, min_deepness=None):
    terms_to_keep = set()

    t_id2level = onto.get_t_id2level()

    def keep(term):
        t_id = term.get_id()
        if t_id in terms_to_keep:
            return
        terms_to_keep.add(t_id)
        for parent_id in term.get_parent_ids():
            if not min_deepness or max(t_id2level[parent_id]) >= min_deepness:
                # add2map(level2term, max(onto.getLevel(p_term)), p_term.getName())
                p_term = onto.get_term(parent_id)
                keep(p_term)
        for (subj, rel, obj) in onto.get_term_relationships(t_id, None, 0):
            if relationships and not (rel in relationships):
                continue
            keep(onto.get_term(subj))
            keep(onto.get_term(obj))

    for t in terms_collection:
        keep(t)

    for term in (t for t in onto.get_all_terms() if t.get_id() not in terms_to_keep):
        onto.remove_term(term, True)

    onto.filter_relationships(relationships)


def save(onto, path):
    processed = set()
    with open(path, 'w') as f:
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
                f.write('synonym: "{0}" RELATED [ChEBI:]\n'.format(syn))
            for parent in term.get_parent_ids():
                f.write("is_a: {0}\n".format(parent))
            for (subj, rel, obj) in onto.get_term_relationships(id_, None, 1):
                f.write("relationship: {1} {0}\n".format(obj, rel))
            f.write("\n")
        for rel in onto.get_relationships():
            f.write("[Typedef]\n")
            f.write("id: {0}\n".format(rel))
            f.write("name: {0}\n".format(rel.replace("_", " ")))
            f.write("is_cyclic: false\n")
            f.write("is_transitive: false\n")


def parse(obo_file, relationships=None):
    if not obo_file or obo_file.find(".obo") == -1 or not os.path.exists(obo_file):
        return None
    ontology = Ontology()
    term = None
    child2parents = {}
    with open(obo_file, 'r') as obo:
        parents = set()
        for line in obo:
            line = line.replace("\n", '')
            if line.find("[Term]") != -1:
                if term:
                    ontology.add_term(term)
                    if parents:
                        child2parents[term.get_id()] = parents
                    parents = set()
                term = Term()
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
                term.set_name(value)
            elif prefix == "is_a":
                parent = ontology.get_term(value)
                if parent:
                    parent.children.add(term)
                else:
                    parents.add(value.lower())
                term.parents.add(value.lower())
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
                end = value.find('"', start + 1)
                if end == -1:
                    continue
                value = value[start + 1:end]
                term.add_synonym(value)
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
                term.add_alt_id(value.strip().replace(" ", "."))
                if value.find("KEGG COMPOUND:") != -1:
                    term.add_kegg(value.replace("KEGG COMPOUND:", '').strip())

    if term:
        ontology.add_term(term)
    for child, parents in child2parents.iteritems():
        child = ontology.get_term(child)
        for parent in parents:
            ontology.get_term(parent).children.add(child)
    return ontology


class Term:
    def __init__(self, t_id=None, name=None, parents=None, children=None):
        self.id = t_id.lower() if t_id else None
        self.altIds = set()
        self.name = name
        self.parents = set(parents) if parents else set()
        self.children = set(children) if children else set()
        self.synonyms = set()
        self.kegg = set()

    def get_id(self):
        return str(self.id)

    def add_alt_id(self, t_id):
        self.altIds.add(t_id.lower())

    def get_all_ids(self):
        return self.altIds | {self.id}

    def set_id(self, t_id):
        self.id = t_id.lower()

    def get_name(self):
        return str(self.name) if self.name else None

    def set_name(self, name):
        self.name = name

    def add_synonym(self, synonym):
        self.synonyms.add(synonym)

    def get_synonyms(self):
        return set(self.synonyms)

    def get_parent_ids(self):
        return set(self.parents)

    def add_parent(self, parent_id):
        self.parents.add(parent_id)

    def add_child(self, term):
        self.children.add(term)

    def get_descendants(self, direct=True):
        result = set(self.children)
        if direct:
            return result
        for child in self.children:
            result |= child.get_descendants(direct)
        return result

    def add_kegg(self, kegg):
        if kegg:
            self.kegg.add(kegg.lower().strip())

    def get_kegg_ids(self):
        return set(self.kegg)

    def __str__(self):
        return "{0}({1})".format(self.get_id(), self.get_name())


def get_name_bis(name):
    return ''.join(e for e in name if e.isalnum())


def is_a(child, parent):
    return child in parent.get_descendants(False)


class Ontology:
    def __init__(self):
        self.roots = set()
        self.id2term = {}
        self.alt_id2term = {}
        self.name2term_ids = defaultdict(set)
        self.rel_map = defaultdict(set)
        self.kegg2term_id = {}

    def get_all_terms(self):
        return self.id2term.values()

    def get_all_term_ids(self):
        return set(self.id2term.iterkeys())

    def add_relationship(self, subj, rel, obj):
        relationship = (subj, rel, obj)
        self.rel_map[subj].add(relationship)
        self.rel_map[obj].add(relationship)
        self.rel_map[rel].add(relationship)

    # role: 0 for any, 1 for subj, 2 for obj
    def get_term_relationships(self, term_id, rel=None, role=0):
        if not (term_id in self.rel_map):
            return iter(())
        relationships = set(self.rel_map[term_id])
        if rel:
            relationships = ((subj, r, obj) for (subj, r, obj) in relationships if rel == r)
        if 1 == role:
            relationships = ((subj, r, obj) for (subj, r, obj) in relationships if term_id == subj)
        if 2 == role:
            relationships = ((subj, r, obj) for (subj, r, obj) in relationships if term_id == obj)
        return relationships

    def get_relationship_participants(self, rel):
        return set(self.rel_map[rel])

    def get_relationships(self):
        result = set()
        for rel_set in self.rel_map.itervalues():
            result |= {rel for (subj, rel, obj) in rel_set}
        return result

    def get_t_id_by_kegg(self, kegg):
        if kegg:
            kegg = kegg.lower().strip()
        if kegg in self.kegg2term_id:
            return self.kegg2term_id[kegg]
        return None

    def add_term(self, term):
        if not term:
            return
        t_id = term.get_id()
        self.id2term[t_id] = term
        for alt_id in term.get_all_ids():
            alt_id = alt_id
            self.alt_id2term[alt_id] = term
        names = set(term.get_synonyms())
        names.add(term.get_name())
        for name in names:
            name = name.lower().strip()
            name_bis = get_name_bis(name)
            self.name2term_ids[name].add(t_id)
            self.name2term_ids[name_bis].add(t_id)
        if not term.get_parent_ids():
            self.roots.add(term)
        for child in term.get_descendants():
            child.parents.add(t_id)
            self.roots -= {child}
        for kegg in term.get_kegg_ids():
            if kegg:
                self.kegg2term_id[kegg.lower().strip()] = t_id

    def filter_relationships(self, rel_to_keep):
        to_remove = set()
        for rel_set in self.rel_map.itervalues():
            for (subj, rel, obj) in rel_set:
                if rel not in rel_to_keep:
                    to_remove.add((subj, rel, obj))
        for (subj, rel, obj) in to_remove:
            if subj in self.rel_map:
                self.rel_map[subj] -= {(subj, rel, obj)}
                if not self.rel_map[subj]:
                    del self.rel_map[subj]
            if obj in self.rel_map:
                self.rel_map[obj] -= {(subj, rel, obj)}
                if not self.rel_map[obj]:
                    del self.rel_map[obj]

    def remove_term(self, term, brutally=False):
        if not term:
            return
        t_id = term.get_id()
        if t_id in self.id2term:
            del self.id2term[t_id]
        for alt_id in term.get_all_ids():
            if alt_id in self.alt_id2term:
                del self.alt_id2term[alt_id]
        names = set(term.get_synonyms())
        names.add(term.get_name())
        for name in names:
            name = name.lower()
            name_bis = get_name_bis(name)
            for it in [name, name_bis]:
                if it in self.name2term_ids:
                    self.name2term_ids[it] -= {t_id}
                    if not self.name2term_ids[it]:
                        del self.name2term_ids[it]
        for kegg in term.get_kegg_ids():
            if kegg in self.kegg2term_id:
                del self.kegg2term_id[kegg]
        parents = term.get_parent_ids()
        if not parents:
            self.roots -= {term}
        children = term.get_descendants()
        for par_id in parents:
            par = self.get_term(par_id)
            if par:
                par.children -= {term}
                if not brutally:
                    par.children |= children
        for child in children:
            child.parents -= term.get_all_ids()
            if not brutally:
                child.parents |= parents
            if not child.parents:
                self.roots.add(child)
        for (subj, rel, obj) in self.get_term_relationships(t_id):
            if t_id == subj and t_id != obj:
                remove_from_map(self.rel_map, obj, (subj, rel, obj))
            elif t_id == obj:
                remove_from_map(self.rel_map, subj, (subj, rel, obj))
        if t_id in self.rel_map:
            del self.rel_map[t_id]

    def get_term(self, term_id):
        if not term_id:
            return None
        term_id = term_id.lower()
        if term_id in self.id2term:
            return self.id2term[term_id]
        if term_id in self.alt_id2term:
            return self.alt_id2term[term_id]
        return None

    def part_of(self, part_id, whole_ids):
        term = self.get_term(part_id)
        if not term:
            return None
        part_of = lambda t: {obj for (subj, r, obj) in self.get_term_relationships(t.get_id(), "part_of", role=1)}
        whole_ids = {t_id.lower().strip() for t_id in whole_ids}
        result = whole_ids & part_of(term)
        if result:
            return result
        term_set = {term}
        result = set()
        part = self.get_term(part_id)
        while term_set:
            items = set()
            for term in term_set:
                parents = term.get_parent_ids()
                wholes = part_of(term)
                candidates = parents | wholes
                for it in candidates:
                    candidate = self.get_term(it)
                    if (it in whole_ids) and not is_a(part, candidate):
                        result.add(it)
                        continue
                    result |= whole_ids & part_of(candidate)
                    result |= {t_id for t_id in set(whole_ids) & candidate.get_parent_ids() if
                               not is_a(part, self.get_term(t_id))}
                    items.add(candidate)
            term_set = items
        return result

    def get_ancestors(self, term, direct=True, rel=None, checked=None):
        if not checked:
            checked = set()
        parents = term.get_parent_ids() if not rel else {obj for (subj, rel, obj) in
                                                         self.get_term_relationships(term.get_id(), rel, 1)}
        direct_parents = set(map(lambda t_id: self.get_term(t_id), parents))
        if direct:
            return direct_parents
        result = set(direct_parents)
        checked.add(term)
        for parent in direct_parents:
            if not (parent in checked):
                result |= self.get_ancestors(parent, direct, rel, checked)
        return result

    def get_t_id2level(self):
        t_id2level = defaultdict(set)

        def _update_level(t, cur_level):
            t_id2level[t.get_id()].add(cur_level)
            for it in t.get_descendants(True):
                _update_level(it, cur_level + 1)

        for t in self.roots:
            _update_level(t, 0)

        return t_id2level

    def get_level(self, term):
        parents = self.get_ancestors(term)
        if not parents:
            return [0]
        level = set()
        for p in parents:
            level |= set(self.get_level(p))
        return [1 + i for i in level]

    def get_equivalents(self, term, rel=None, direction=0, relationships=None):
        term_id = term.get_id()
        equals = set()
        for (subj, r, obj) in self.get_term_relationships(term_id, rel, direction):
            if not relationships or r in relationships:
                equals.add(obj if subj == term_id else subj)
        return {self.get_term(t_id) for t_id in equals}

    def get_sub_tree(self, t, relationships=None):
        return self.get_generalized_descendants(t, False, set(), relationships) \
               | self.get_equivalents(t, None, 0, relationships) | {t}

    def get_generalized_descendants(self, term, direct=True, checked=None, relationships=None):
        if checked is None:
            checked = set()
        term2eqs = lambda t: {t} | self.get_equivalents(t, None, 0, relationships)
        cup = lambda s1, s2: s1 | s2
        terms = term2eqs(term)
        direct_children = \
            reduce(cup, (reduce(cup, (term2eqs(t) for t in it.get_descendants(True)), set()) for it in terms), set())
        if direct:
            return direct_children
        checked |= terms
        return reduce(cup, (self.get_generalized_descendants(t, direct, checked, relationships)
                            for t in direct_children - checked), direct_children)

    def get_generalized_ancestors(self, term, direct=True, checked=None, relationships=None, depth=None):
        if depth is not None and depth <= 0:
            return set()
        if not checked:
            checked = set()
        terms = {term} | self.get_equivalents(term, None, 0, relationships)
        direct_parents = set()
        for it in terms:
            parents = {self.get_term(t_id) for t_id in it.get_parent_ids()}
            direct_parents |= parents
            for par in parents:
                direct_parents |= self.get_equivalents(par, None, 0, relationships)
        if direct or 1 == depth:
            return direct_parents
        checked |= terms
        result = set(direct_parents)
        for parent in direct_parents - checked:
            result |= self.get_generalized_ancestors(parent, direct, checked, relationships,
                                                     depth - 1 if depth is not None else depth)
        return result

    def get_generalized_ancestors_of_level(self, term, checked=None, relationships=None, depth=None):
        if depth is None:
            return self.get_generalized_ancestors(term, False, checked, relationships, depth)
        if depth <= 0:
            return set()
        if not checked:
            checked = set()
        terms = {term} | self.get_equivalents(term, None, 0, relationships)
        direct_parents = set()
        for it in terms:
            parents = {self.get_term(t_id) for t_id in it.get_parent_ids()}
            direct_parents |= parents
            for par in parents:
                direct_parents |= self.get_equivalents(par, None, 0, relationships)
        if not direct_parents:
            return terms
        if 1 == depth:
            return direct_parents
        checked |= terms
        result = set()
        for parent in direct_parents - checked:
            result |= self.get_generalized_ancestors_of_level(parent, checked, relationships, depth - 1)
        return result

    def get_roots(self):
        return set(self.roots)

    def get_ids_by_name(self, name):
        name = name.lower()
        name_bis = get_name_bis(name)
        return (set(self.name2term_ids[name]) if name in self.name2term_ids else set()) \
               | (set(self.name2term_ids[name_bis]) if name_bis in self.name2term_ids else set())

    def common_points(self, terms, depth=None, relationships=None):
        if not terms or depth is not None and depth <= 0:
            return None
        terms = set(terms)
        first = terms.pop()
        common = self.get_generalized_ancestors(term=first, direct=False, checked=set(), depth=depth, relationships=relationships) | \
                 self.get_equivalents(first, relationships=relationships) | {
            first}
        for t in terms:
            common &= self.get_generalized_ancestors(t, False, set(), depth=depth, relationships=relationships) | \
                      self.get_equivalents(t, relationships=relationships) | {t}
        result = set(common)
        return [it for it in common
                if not self.get_generalized_descendants(it, False, set(), relationships=relationships) & result]

    def remove_relationships(self, relationships, brutally=False):
        for (subj_id, r, o_id) in relationships:
            if "is_a" == r:
                subj, obj = self.get_term(subj_id), self.get_term(o_id)
                if not subj or not obj:
                    continue
                subj.parents -= obj.getAllIds()
                obj.children -= {subj}
                if not brutally:
                    subj.parents |= obj.parents
                    for par in obj.parents:
                        self.get_term(par).children.add(subj)
                if not subj.parents:
                    self.roots.add(subj)
            else:
                remove_from_map(self.rel_map, subj_id, (subj_id, r, o_id))
                remove_from_map(self.rel_map, o_id, (subj_id, r, o_id))

    def trim(self, root_ids, relationships=None):
        for r_id in root_ids:
            r = self.get_term(r_id)
            ancestors = self.get_generalized_ancestors(r, relationships=relationships)
            if not {a.id for a in ancestors} & root_ids:
                for it in ancestors:
                    self.remove_term(it, True)
