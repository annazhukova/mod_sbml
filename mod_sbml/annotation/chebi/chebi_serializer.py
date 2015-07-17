import os

from mod_sbml.annotation.chebi.chebi_annotator import EQUIVALENT_RELATIONSHIPS, MOLECULAR_ENTITY, HAS_ROLE_RELATIONSHIP, \
    add_equivalent_chebi_ids, get_cofactor_ids
from mod_sbml.onto import parse, save_simple, save
from mod_sbml import annotation

__author__ = 'anna'

CHEBI = os.path.join(os.path.dirname(os.path.abspath(annotation.__file__)), '..', 'data', 'chebi.obo')
SUB_CHEBI = os.path.join(os.path.dirname(os.path.abspath(annotation.__file__)), '..', 'data', 'subchebi.obo')
SUB_CHEBI_TXT = os.path.join(os.path.dirname(os.path.abspath(annotation.__file__)), '..', 'data', 'subchebi.txt')

COFACTORS_FILE = os.path.join(os.path.dirname(os.path.abspath(annotation.__file__)), '..', 'data', 'ub_cofactors.txt')
COMMON_TERMS_FILE = os.path.join(os.path.dirname(os.path.abspath(annotation.__file__)), '..', 'data', 'ub_common.txt')


# most common ones, like water, H+, oxygen, NAD, etc.
COMMON_UB_IDS = {'chebi:17808', 'chebi:37568', 'chebi:24636', 'chebi:15422', 'chebi:15377', 'chebi:15378',
                 'chebi:15379', 'chebi:37563', 'chebi:37565', 'chebi:17877', 'chebi:18245', 'chebi:16311',
                 'chebi:16192', 'chebi:15846', 'chebi:30089', 'chebi:16497', 'chebi:25805', 'chebi:26020',
                 'chebi:13390', 'chebi:13392', 'chebi:43474', 'chebi:57783', 'chebi:17625', 'chebi:18359',
                 'chebi:25524', 'chebi:58115', 'chebi:16695', 'chebi:15346', 'chebi:17552', 'chebi:58107',
                 'chebi:17330', 'chebi:57692', 'chebi:17013', 'chebi:61404', 'chebi:30616', 'chebi:61402',
                 'chebi:57299', 'chebi:58307', 'chebi:58223', 'chebi:35780', 'chebi:17659', 'chebi:28862',
                 'chebi:58342', 'chebi:57287', 'chebi:17093', 'chebi:16027', 'chebi:58189', 'chebi:33813',
                 'chebi:29888', 'chebi:16908', 'chebi:57540', 'chebi:456216', 'chebi:456215', 'chebi:57945',
                 'chebi:28971', 'chebi:16474', 'chebi:28850', 'chebi:18421', 'chebi:16761', 'chebi:17361',
                 'chebi:26689', 'chebi:57288', 'chebi:16284', 'chebi:16234', 'chebi:16174', 'chebi:82680',
                 'chebi:18075', 'chebi:16171', 'chebi:18077', 'chebi:15713', 'chebi:16039', 'chebi:16238',
                 'chebi:36080', 'chebi:29375', 'chebi:73342', 'chebi:15366', 'chebi:16526', 'chebi:17984',
                 'chebi:17544', 'chebi:28846', 'chebi:29033', 'chebi:17677', 'chebi:15996', 'chebi:18009',
                 'chebi:58210', 'chebi:61429', 'chebi:58280', 'chebi:18361', 'chebi:15918', 'chebi:246422',
                 'chebi:17239', 'chebi:73627', 'chebi:13389', 'chebi:16240', 'chebi:58245', 'chebi:29325',
                 'chebi:13534', 'chebi:15351', 'chebi:26078'}


def get_chebi():
    return SUB_CHEBI_TXT


def get_closest_descendants_with_formulas(onto, term):
    if term.get_formulas():
        return {term}
    for t in onto.get_equivalents(term, relationships=EQUIVALENT_RELATIONSHIPS):
        if t.get_formulas():
            return {term}
    return reduce(lambda s1, s2: s1 | s2,
                  (get_closest_descendants_with_formulas(onto, t)
                   for t in onto.get_generalized_descendants(term, True)),
                  set())


def serialize_sub_chebi():
    onto = parse(CHEBI, EQUIVALENT_RELATIONSHIPS)
    terms_to_keep = reduce(lambda s1, s2: s1 | s2,
                           (onto.get_sub_tree(t, relationships=EQUIVALENT_RELATIONSHIPS) for t in
                               get_closest_descendants_with_formulas(onto, onto.get_term(MOLECULAR_ENTITY))), set())
    for term in (t for t in onto.get_all_terms() if t not in terms_to_keep):
        onto.remove_term(term, True)
    save_simple(onto, SUB_CHEBI_TXT)
    save(onto, SUB_CHEBI)


def serialize_ubiquitous_terms():
    chebi = parse(CHEBI, EQUIVALENT_RELATIONSHIPS | {HAS_ROLE_RELATIONSHIP})
    ub_ch_ids = add_equivalent_chebi_ids(chebi, COMMON_UB_IDS)
    with open(COMMON_TERMS_FILE, 'w+') as f:
        f.write('\t'.join(ub_ch_ids))

    ub_ch_ids = get_cofactor_ids(chebi)
    with open(COFACTORS_FILE, 'w+') as f:
        f.write('\t'.join(ub_ch_ids))


# serialize_sub_chebi()
# serialize_ubiquitous_terms()
