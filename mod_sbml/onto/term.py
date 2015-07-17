__author__ = 'anna'


class Term:
    def __init__(self, onto, t_id=None, name=None, parent_ids=None, s=None):
        self.onto = onto
        self.id = t_id.lower() if t_id else None
        self.altIds = set()
        self.name = name
        self.synonyms = set()
        self.kegg = set()
        self.parent_ids = set()
        self.formulas = set()
        if parent_ids:
            for p_id in parent_ids:
                self.add_parent(p_id)
        if s:
            self.__from_string(s)

    def get_id(self):
        return str(self.id)

    def add_alt_id(self, t_id):
        self.altIds.add(t_id.lower())

    def get_all_ids(self):
        return self.altIds | {self.id}

    def set_id(self, t_id):
        if t_id:
            self.id = t_id.lower()

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

    def add_formula(self, formula):
        formula = formula.strip()
        if formula and formula != '.':
            self.formulas.add(formula)

    def get_synonyms(self):
        return set(self.synonyms)

    def get_formulas(self):
        return set(self.formulas)

    def get_parent_ids(self):
        return set(self.parent_ids)

    def add_parent(self, parent_id):
        parent_id = parent_id.strip()
        if parent_id:
            parent_id = parent_id.lower()
            self.parent_ids.add(parent_id)
            self.onto.parent2children[parent_id] |= {self.id}

    def add_kegg(self, kegg):
        if kegg:
            self.kegg.add(kegg.lower().strip())

    def get_kegg_ids(self):
        return set(self.kegg)

    def __from_string(self, s):
        if not s:
            return
        s = s.splitlines()[0]
        t_id, alt_ids, name, synonyms, parent_ids, formulas = s.split('\t\t')
        self.set_id(t_id)
        self.set_name(name)
        for synonym in synonyms.split('\t'):
            self.add_synonym(synonym)
        for alt_id in alt_ids.split('\t'):
            self.add_alt_id(alt_id)
        for parent_id in parent_ids.split('\t'):
            self.add_parent(parent_id)
        for formula in formulas.split('\t'):
            self.add_formula(formula)

    def __str__(self):
        return "%s\t\t%s\t\t%s\t\t%s\t\t%s\t\t%s\n" % (self.id, '\t'.join(self.altIds),
                                                       self.name, '\t'.join(self.synonyms),
                                                       '\t'.join(self.parent_ids),
                                                       '\t'.join(self.formulas))

    def __eq__(self, other):
        if not other or not isinstance(other, Term):
            return False
        return other.id == self.id

    def __ne__(self, other):
        return not self.__eq__(other)

    def __hash__(self):
        return hash(self.id)

