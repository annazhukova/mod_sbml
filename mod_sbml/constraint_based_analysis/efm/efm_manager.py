import logging
import os

import libsbml

from mod_sbml.constraint_based_analysis.efm.EFM import get_int_size, EFM
from mod_sbml.gibbs.reaction_boundary_manager import get_bounds
from mod_sbml.sbml.sbml_manager import get_products, get_reactants

PRECISION = 6

__author__ = 'anna'


def compute_efms(sbml, directory, em_number, r_id, rev, tree_efm_path, r_id2rev_2threshold=None, threshold=0.0,
                 r_ids=None, rewrite=True):
    """
    Computes elementary flux modes (EFMs) in a given SBML (see http://sbml.org) model,
    that contain a reaction of interest
    (using TreeEFM software [Pey et al. 2014, PMID: 25380956]).

    :param sbml: string, path to the SBML file with the model.
    :param directory: string, directory where to store the results, such as stoichiometric matrix, EFMs.
    :param tree_efm_path: string,path to the executable of TreeEFM software [Pey et al. 2014, PMID: 25380956].
    :param r_id:string, id of the reaction of interest.
    :param rev: boolean, if the reaction of interest should be considered in the opposite direction.
    :param em_number: int, number of EFMs to compute with TreeEFM software [Pey et al. 2014, PMID: 25380956].
    :param r_id2rev_2threshold: set of strings in the form {(r_id_0, reversed_0), (r_id_1, reversed_1), ...},
    if specified, only EFMs that contain all (if reaction_op == (REACTION_OPERATION_AND))
    or any (if reaction_op == (REACTION_OPERATION_OR) of the reaction ids in specified directions
    from this set will be returned.
    :return:efms: list of EFMs

    :raise ValueError: if the reaction of interest was not found in the model.
    """
    if rewrite and not os.path.exists(tree_efm_path):
        raise ValueError("TreeEFM runner is not found at %s" % tree_efm_path)
    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()
    # Compute stoichiometric matrix
    st_matrix_file = os.path.join(directory, "st_matrix.txt")
    s_id2i, r_id2i, rev_r_id2i = stoichiometric_matrix(model, st_matrix_file)
    logging.info("stoichiometric matrix saved to %s" % st_matrix_file)
    # Figure out in which reaction we are interested in
    if (rev and r_id not in rev_r_id2i) or (not rev and r_id not in r_id2i):
        raise ValueError("R%seaction with id %s is not found in the model" % ('eversed r' if rev else '', r_id))
    i = rev_r_id2i[r_id] if rev else r_id2i[r_id]
    # Compute EFMs using TreeEFM software (Pey et al. 2014, PMID: 25380956)
    efm_file = os.path.join(directory, "FV-EM.dat")
    efm_file_txt = "%s.txt" % efm_file
    if rewrite or not os.path.exists(efm_file_txt):
        os.system("%s -r %d -i %s -l EM -e %d -o %s -z %d" % (tree_efm_path, i, st_matrix_file, em_number, directory,
                                                              threshold))
        os.system("%s -b %s" % (tree_efm_path, efm_file))
        logging.info("elementary modes saved to %s" % efm_file_txt)
    # Filter EFMs so that only those that don't include the reaction in opposite directions are left.
    # If r_id2rev are specified, filter EFMs to leave only those that include these reactions in these directions.
    em_file_filtered = os.path.join(directory, "FV-EM_filtered.dat.txt")
    efms, all_r_ids, rev_r_ids = filter_efms(efm_file_txt, r_id2i, rev_r_id2i, em_file_filtered, r_id2rev_2threshold,
                                             zero_threshold=threshold, r_ids=r_ids)
    logging.info(
        "%d elementary modes corresponding to reactions of interest saved to %s" % (len(efms), em_file_filtered))
    efms = sorted(efms, key=lambda efm: len(efm))
    return efms, all_r_ids, rev_r_ids


def filter_efms(in_path, r_id2i, rev_r_id2i, out_path, r_id2rev_2threshold, zero_threshold=0.0,
                r_ids=None):
    i2r_id = {i: (r_id, False) for (r_id, i) in r_id2i.iteritems()}
    i2r_id.update({i: (r_id, True) for (r_id, i) in rev_r_id2i.iteritems()})
    all_r_ids = r_ids if r_ids else sorted(set(r_id2i.iterkeys()) | set(rev_r_id2i.iterkeys()))
    rev_r_ids = set(rev_r_id2i.iterkeys())
    int_size = get_int_size()
    processed = set()
    efms = []
    round_value = lambda v: round(float(v), PRECISION)
    rejected_bad, rejected_different = 0, 0
    with open(out_path, 'w+') as out_f:
        with open(in_path, 'r') as in_f:
            for line in in_f:
                values = line.replace("\n", "").strip().split(" ")
                r_id2coefficient = {}
                bad_em = False
                for (v, i) in zip(values, xrange(1, len(values) + 1)):
                    v = round_value(v)
                    if not v or abs(v) <= zero_threshold:
                        continue
                    r_id, rev = i2r_id[i]
                    if rev:
                        v *= -1
                    # The same reaction participates in different directions
                    # => don't want such an EFM
                    if r_id in r_id2coefficient:
                        bad_em = True
                        break
                    r_id2coefficient[r_id] = v

                # The same reaction participates in different directions
                # => don't want such an EFM
                if bad_em:
                    rejected_bad += 1
                    continue

                if r_id2rev_2threshold:
                    for (r_id2rev, present_reaction_threshold) in r_id2rev_2threshold:
                        present_r_ids = set()
                        for (r_id, rev) in r_id2rev.iteritems():
                            if r_id in r_id2coefficient and (rev is None or rev == (r_id2coefficient[r_id] < 0)):
                                present_r_ids.add(r_id)
                        if len(present_r_ids) * 1.0 / len(r_id2rev) < present_reaction_threshold:
                            rejected_different += 1
                            bad_em = True
                            break
                    if bad_em:
                        continue

                out_f.write(line)

                if r_ids:
                    r_id2coefficient = {r_id: val for (r_id, val) in r_id2coefficient.iteritems() if r_id in r_ids}
                efm = EFM(r_id2coeff=r_id2coefficient, r_ids=all_r_ids, rev_r_ids=rev_r_ids, int_size=int_size)

                if not r_ids or efm.binary_efm not in processed:
                    efms.append(efm)
                    if r_ids:
                        processed.add(efm.binary_efm)
    if rejected_different:
        logging.info('Rejected %d EFMs as not all of the reactions of interest were present in them.'
                     % rejected_different)
    if rejected_bad:
        logging.info('Rejected %d EFMs as they contained reversible reactions in both directions' % rejected_bad)
    return efms, all_r_ids, rev_r_ids


def stoichiometric_matrix(model, path):
    """
    Extracts the model's stoichiometric matrix in a format compatible with TreeEFM [Pey et al., 2014],
    i.e. one cell per line as follows: row,column,value.
    :param model: libsbml model
    :param path: path to file where to save the matrix
    :return m_id2i, r_id2i, rev_r_id2i: a tuple of three dictionaries: species id to its index,
    reactions id to its indices; reversible reaction id to its index (for the opposite direction).
    """
    internal_s_ids = [s.id for s in model.getListOfSpecies() if not s.getBoundaryCondition()]
    m_id2i = dict(zip(internal_s_ids, xrange(1, len(internal_s_ids) + 1)))
    r_id2i, rev_r_id2i = {}, {}

    def add_reaction_data(file, reaction_number, reaction, rev=False):
        for (m_id, st) in get_reactants(reaction, True):
            if m_id in m_id2i:
                file.write("%d,%d,%d\n" % (m_id2i[m_id], reaction_number, st if rev else -st))
        for (m_id, st) in get_products(reaction, True):
            if m_id in m_id2i:
                file.write("%d,%d,%d\n" % (m_id2i[m_id], reaction_number, -st if rev else st))

    i = 1
    with open(path, 'w+') as f:
        for r in sorted(model.getListOfReactions(), key=lambda r: r.id):
            l_b, u_b = get_bounds(r)
            if u_b > 0:
                add_reaction_data(f, i, r)
                r_id2i[r.id] = i
                i += 1
            if l_b < 0:
                add_reaction_data(f, i, r, rev=True)
                rev_r_id2i[r.id] = i
                i += 1
    return m_id2i, r_id2i, rev_r_id2i