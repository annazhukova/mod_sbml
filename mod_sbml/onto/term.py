from collections import defaultdict

__author__ = 'anna'


FORMULA = 'formula'
KEGG = 'kegg.compound'


class Term:
    def __init__(self, onto, t_id=None, name=None, parent_ids=None, s=None):
        self.onto = onto
        self.id = normalize(t_id)
        self.altIds = set()
        self.name = name
        self.synonyms = set()
        self.parent_ids = set()
        self.xrefs = defaultdict(set)
        if parent_ids:
            for p_id in parent_ids:
                self.add_parent(p_id)
        if s:
            self.__from_string(s)

    def get_id(self):
        return str(self.id)

    def add_alt_id(self, t_id):
        t_id = normalize(t_id)
        if t_id:
            self.altIds.add(t_id)

    def get_all_ids(self):
        return self.altIds | {self.id}

    def set_id(self, t_id):
        t_id = normalize(t_id)
        if t_id:
            self.id = t_id

    def get_name(self):
        return str(self.name) if self.name else None

    def set_name(self, name):
        name = name.strip()
        if name:
            self.name = name

    def add_synonym(self, synonym):
        synonym = synonym.strip()
        if synonym:
            self.synonyms.add(synonym)

    def get_synonyms(self):
        return set(self.synonyms)

    def get_formulas(self):
        return self.get_xrefs(FORMULA)

    def get_parent_ids(self):
        return set(self.parent_ids)

    def add_parent(self, parent_id):
        parent_id = normalize(parent_id)
        if parent_id:
            self.parent_ids.add(parent_id)
            self.onto.parent2children[parent_id] |= {self.id}

    def get_kegg_ids(self):
        return self.get_xrefs(KEGG)

    def add_xref(self, db, value):
        db = normalize(db).replace(' ', '.')
        value = value.strip()
        if db and value:
            self.xrefs[db].add(value)

    def get_xrefs(self, db):
        db = normalize(db).replace(' ', '.')
        return set(self.xrefs[db]) if db else set()

    def get_dbs(self):
        return set(self.xrefs.keys())

    def __from_string(self, s):
        if not s:
            return
        s = s.splitlines()[0]
        t_id, alt_ids, name, synonyms, parent_ids, xrefs = s.split('\t\t')
        self.set_id(t_id)
        self.set_name(name)
        for synonym in synonyms.split('\t'):
            self.add_synonym(synonym)
        for alt_id in alt_ids.split('\t'):
            self.add_alt_id(alt_id)
        for parent_id in parent_ids.split('\t'):
            self.add_parent(parent_id)
        for db_values in xrefs.split(';\t;'):
            values = db_values.split('\t')
            if len(values) > 1:
                db = values[0]
                for value in values[1:]:
                    self.add_xref(db, value)

    def __str__(self):
        return "%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\n" % (self.id, '\t'.join(self.altIds),
                                                       self.name, '\t'.join(self.synonyms),
                                                       '\t'.join(self.parent_ids),
                                                       ';\t;'.join(
                                                           ('%s\t%s' % (db, '\t'.join(values))
                                                            for (db, values) in self.xrefs.items())))

    def __eq__(self, other):
        if not other or not isinstance(other, Term):
            return False
        return other.id == self.id

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.id)


def normalize(s):
    if s:
        return s.strip().lower()
    return s

