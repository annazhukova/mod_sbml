from collections import defaultdict

from gibbs.reaction_boundary_manager import get_bounds
from sbml.sbml_manager import get_reactants, get_products


__author__ = 'anna'


def stoichiometric_matrix(model, path):
    """
    Extracts the model's stoichiometric matrix in a format compatible with TreeEFM [Pey et al., 2014],
    i.e. one cell per line as follows: row,column,value.
    :param model: libsbml model
    :param path: path to file where to save the matrix
    :return s_id2i, r_id2i: a tuple of dictionaries: species id to its index,
    reactions id to its indices: (i, i+1) for reversible reactions, (i, 0) for irreversible ones.
    """
    internal_s_ids = [s.id for s in model.getListOfSpecies() if not s.getBoundaryCondition()]
    s_id2i = dict(zip(internal_s_ids, xrange(1, len(internal_s_ids) + 1)))
    r_id2i = {}
    i = 1
    with open(path, 'w+') as f:
        for r in model.getListOfReactions():
            l, u = 0, 0
            l_b, u_b = get_bounds(r)
            if u_b > 0:
                l = i
                i += 1
                add_reaction_data(f, l, r, s_id2i)
            if l_b < 0:
                u = i
                i += 1
                add_reaction_data(f, u, r, s_id2i, reversed=True)
            r_id2i[r.id] = (l, u)
    return s_id2i, r_id2i


def add_reaction_data(file, reaction_number, reaction, m_id2i, reversed=False):
    for (m_id, st) in get_reactants(reaction, True):
        if m_id in m_id2i:
            file.write("%d,%d,%d\n" % (m_id2i[m_id], reaction_number, st if reversed else -st))
    for (m_id, st) in get_products(reaction, True):
        if m_id in m_id2i:
            file.write("%d,%d,%d\n" % (m_id2i[m_id], reaction_number, -st if reversed else st))


# def ems2binary(in_path, r_id2i):
#     r_ids = []
#     for r_id, (i, j) in sorted(r_id2i.iteritems(), key=lambda (r_id, ij): max(ij)):
#         if i:
#             r_ids.append(r_id)
#         if j:
#             r_ids.append("%s_rev" % r_id)
#
#     binary_ems = []
#     with open(in_path, 'r') as in_f:
#         for line in in_f:
#             i = 1
#             em = []
#             line = line.replace("\n", "").strip()
#             values = line.split(" ")
#             for v in values:
#                 fv = round(float(v), 6)
#                 em.append(1 if fv > 0 else (0 if abs(fv) == 0 else -1))
#                 i += 1
#             binary_ems.append(tuple(em))
#     return binary_ems, r_ids


def ems2binary(efms, r_ids):
    binary_ems = []
    for r_id2value in efms:
        binary_ems.append([(1 if (r_id2value[r_id] > 0) else -1) if r_id in r_id2value else 0 for r_id in r_ids])
    return binary_ems


def format_ems(in_path, r_id2i, threshold=0, r_ids_to_keep=None):
    r_ids = []
    for r_id, (i, j) in sorted(r_id2i.iteritems(), key=lambda (r_id, ij): max(ij)):
        if i:
            r_ids.append((r_id, True))
        if j:
            r_ids.append((r_id, False))
    processed = set()
    ems = []
    with open(in_path, 'r') as in_f:
        for line in in_f:
            em_support = []
            em = defaultdict(list)
            line = line.replace("\n", "").strip()
            values = line.split(" ")
            for v, (r_id, correct_direction) in zip(values, r_ids):
                if r_ids_to_keep and r_id not in r_ids_to_keep:
                    continue
                fv = round(float(v), 6)
                if abs(fv) > threshold:
                    em[r_id] = fv if correct_direction else -fv
                    em_support.append((r_id, correct_direction))
            em_support = tuple(sorted(em_support))
            if em_support not in processed:
                ems.append(em)
                processed.add(em_support)
    return ems


def save_ems_for_acom(ems, react_file, dat_file):
    rs = set()
    r_id_ems = []
    i = 0
    for em in ems:
        r_id_em = [(r_id, '1' if sum(coeffs) > 0 else '0') for (r_id, coeffs) in em if sum(coeffs) != 0]
        rs |= {it[0] for it in r_id_em}
        r_id_ems.append(r_id_em)
        i += 1
        if i > 10:
            break
    rs = sorted(rs)
    with open(react_file, 'w+') as f:
        f.write(",".join(rs))
    with open(dat_file, 'w+') as f:
        for em in r_id_ems:
            f.write("%s\n" % " ".join((em[r] if r in em else '0' for r in rs)))







