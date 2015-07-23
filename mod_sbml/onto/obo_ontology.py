#!/usr/bin/env python
# encoding: utf-8
from collections import defaultdict

from mod_sbml.utils.misc import remove_from_map
from mod_sbml.utils.natsort import natsorted, natcasecmp

PART_OF = "part_of"

__author__ = 'anna'


def normalize(name):
    return ''.join(e for e in name if e.isalnum()).lower()


class Ontology:
    def __init__(self):
        self.roots = set()
        self.id2term = {}
        self.alt_id2term = {}
        self.name2term_ids = defaultdict(set)
        self.rel_map = defaultdict(set)
        self.xref2term_ids = defaultdict(set)
        self.parent2children = defaultdict(set)

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
            name = normalize(name)
            if name:
                self.name2term_ids[name].add(t_id)
        if not term.get_parent_ids():
            self.roots.add(term)
        for db in term.get_dbs():
            for value in term.get_xrefs(db):
                value = value.lower()
                self.xref2term_ids[value].add(t_id)
                self.xref2term_ids['%s:%s' % (db, value)].add(t_id)

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
            name = normalize(name)
            if name and name in self.name2term_ids:
                self.name2term_ids[name] -= {t_id}
                if not self.name2term_ids[name]:
                    del self.name2term_ids[name]
        for db in term.get_dbs():
            for value in term.get_xrefs(db):
                value = value.lower()
                self.xref2term_ids[value] -= {t_id}
                if not self.xref2term_ids[value]:
                    del self.xref2term_ids[value]
                db_value = '%s:%s' % (db, value)
                self.xref2term_ids[db_value] -= {t_id}
                if not self.xref2term_ids[db_value]:
                    del self.xref2term_ids[db_value]
        parents = term.get_parent_ids()
        if not parents:
            self.roots -= {term}
        child_ids = self.get_descendants(t_id)
        for child_id in child_ids:
            child = self.get_term(child_id)
            if not child:
                continue
            child.parent_ids -= term.get_all_ids()
            if not brutally:
                child.parent_ids |= parents
            if not child.parent_ids:
                self.roots.add(child)
        if t_id in self.parent2children:
            del self.parent2children[t_id]
            for par_id in parents:
                self.parent2children[par_id] -= {t_id}
                if not brutally:
                    self.parent2children[par_id] |= child_ids

        for (subj, rel, obj) in self.get_term_relationships(t_id):
            if t_id == subj and t_id != obj:
                remove_from_map(self.rel_map, obj, (subj, rel, obj))
            elif t_id == obj:
                remove_from_map(self.rel_map, subj, (subj, rel, obj))
        if t_id in self.rel_map:
            del self.rel_map[t_id]

    def get_term(self, key, check_only_ids=True):
        """
        Looks for a term corresponding to the given key.
        By default, the key is treated as an id or alternative id.
        If check_only_ids argument is set to False (by default it's True),
        the key is also looked for in term names and xrefs.
        :param key: str, by default the term's id or alternative id.
        If check_only_ids argument is set to False (by default it's True),
        the key is also looked for in term names and xrefs.
        :param check_only_ids: boolean, optional. If set to False (by default it's True),
        the key is also looked for in term names and xrefs.
        :return: term (instance of class mod_sbml.onto.term.Term) corresponding to the given key,
        or None if no such term was found.
        """
        if not key:
            return None
        key = key.lower().strip()
        if key in self.id2term:
            return self.id2term[key]
        if key in self.alt_id2term:
            return self.alt_id2term[key]
        if not check_only_ids:
            if key in self.xref2term_ids:
                for t_id in natsorted(self.xref2term_ids[key], cmp=natcasecmp):
                    if t_id in self.id2term:
                        return self.id2term[t_id]
            key = normalize(key)
            if key in self.name2term_ids:
                for t_id in natsorted(self.name2term_ids[key], cmp=natcasecmp):
                    if t_id in self.id2term:
                        return self.id2term[t_id]
        return None

    def get_descendants(self, term_id, direct=True):
        result = set(self.parent2children[term_id])
        if direct:
            return result
        for child_id in self.parent2children[term_id]:
            result |= self.get_descendants(child_id, direct)
        return result

    def is_a(self, child_id, parent_id):
        return child_id and child_id.lower() in self.get_descendants(parent_id, False)

    def part_of(self, part_id, whole_ids):
        part = self.get_term(part_id)
        if not part:
            return None
        part_of_rel = lambda t_id: {obj for (subj, r, obj) in self.get_term_relationships(t_id, PART_OF, role=1)}
        whole_ids = {t_id.lower().strip() for t_id in whole_ids}
        result = whole_ids & part_of_rel(part_id)
        if result:
            return result
        term_set = {part}
        result = set()
        while term_set:
            items = set()
            for term in term_set:
                parents = term.get_parent_ids()
                wholes = part_of_rel(term.get_id())
                candidates = parents | wholes
                for candidate_id in candidates:
                    if (candidate_id in whole_ids) and not self.is_a(part_id, candidate_id):
                        result.add(candidate_id)
                        continue
                    candidate = self.get_term(candidate_id)
                    result |= whole_ids & part_of_rel(candidate_id)
                    result |= {t_id for t_id in set(whole_ids) & candidate.get_parent_ids() if
                               not self.is_a(part_id, t_id)}
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

        def _update_level(t_id, cur_level):
            t_id2level[t_id].add(cur_level)
            for it in self.get_descendants(t_id, True):
                _update_level(it, cur_level + 1)

        for t in self.roots:
            _update_level(t.get_id(), 0)

        return t_id2level

    def get_level(self, term):
        parents = self.get_ancestors(term)
        if not parents:
            return [0]
        level = set()
        for p in parents:
            level |= set(self.get_level(p))
        return [1 + i for i in level]

    def get_equivalents(self, term, rel=None, direction=0, relationships=None, checked=None):
        term_id = term.get_id()
        if checked is None:
            checked = set()
        checked.add(term_id)
        equals = set()
        for (subj, r, obj) in self.get_term_relationships(term_id, rel, direction):
            if not relationships or r in relationships:
                eq_t_id = obj if subj == term_id else subj
                if eq_t_id in checked:
                    continue
                checked.add(eq_t_id)
                eq_term = self.get_term(eq_t_id)
                equals.add(eq_term)
                equals |= self.get_equivalents(eq_term, rel, direction, relationships, checked)
        return equals

    def get_sub_tree(self, t, relationships=None, depth=None):
        return self.get_generalized_descendants(t, False, set(), relationships, depth=depth) \
               | self.get_equivalents(t, None, 0, relationships) | {t}

    def get_generalized_descendants(self, term, direct=True, checked=None, relationships=None, depth=None):
        if depth is not None and depth <= 0:
            return set()
        if checked is None:
            checked = set()
        term2eqs = lambda t: {t} | self.get_equivalents(t, None, 0, relationships)
        cup = lambda s1, s2: s1 | s2
        terms = term2eqs(term)
        direct_children = \
            reduce(cup,
                   (reduce(cup,
                           (term2eqs(self.get_term(t)) for t in self.get_descendants(it.get_id(), True)),
                           set()) for it in terms),
                   set())
        if direct or 1 == depth:
            return direct_children
        checked |= terms
        return reduce(cup, (self.get_generalized_descendants(t, direct, checked, relationships,
                                                             depth - 1 if depth is not None else depth)
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
