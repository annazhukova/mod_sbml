import logging
import os
import math
import sys
import gmpy

import libsbml

from em.acom import classify
from em.stoichiometry_manager import stoichiometric_matrix, ems2binary

PRECISION = 6

__author__ = 'anna'


def binary2efm(binary_efm, r_ids, reversible_r_ids, int_size=None, coefficients=None, binary=False):
    """
    Converts an EFM from a binary representation into a dictionary
    that maps ids of active reactions to their coefficients: {r_id: coefficient}.
    If binary is set to True, the coefficient values are 1 (for any active reaction in the standard direction)
    or -1 (for reactions that are active in the reversed direction).

    :param binary_efm: an EFM in a binary representation.
    A binary representation of an EFM is a list of integers whose binary representations
    correspond to the reactions that are active in the EFM: if the reaction is active,
    the corresponding bit is set to 1.
    If the total number of reactions in the model is larger that the number of bits in an int, several ints are used.

    Example: For a model containing reactions r1, r2(reversible), r3(reversible), r4, r5,
    a EFM: 3 r1, -2 r2, 1 r3, 1 r5 would be represented as [77], as the binary representation of 77 is '1001101'
    that corresponds to '1 r5, 0 r4, 0 -r3, 1 r3, 1 -r2, 0 r2, 1 r1', where -ri is the reversed version of ri.


    :param r_ids: ordered collection of reaction ids (strings).

    :param reversible_r_ids: set of ids of reversible reactions (strings).

    :param int_size: the maximal number of bits in an int, can be calculated as math.log(sys.maxint) / math.log(2).

    :param coefficients: (optional) a list of non-zero coefficients corresponding to the binary EFM.
    If not present, the coefficients are assumed to take values 1 or -1.
    Example: The non_zero_coefficients for the EFM in the example above are [3, -2, 1, 1].

    :param binary: boolean, if is set to True, the coefficient values in the result will be
    1 (for any active reaction in the standard direction)
    or -1 (for reactions that are active in the reversed direction).
    Otherwise, any float coefficient values will be possible.

    :return:
    """
    if not int_size:
        int_size = get_int_size()
    converted_efm = {}
    binary_efm_iterable = iter(binary_efm)
    cur_efm_part = next(binary_efm_iterable)
    coeff_iterable = iter(coefficients) if coefficients else None

    def process(r_id, cur_efm_part, i, reversed=False):
        if cur_efm_part & i:
            coeff = next(coeff_iterable) if coefficients else (-1 if reversed else 1)
            if binary:
                converted_efm[r_id] = 1 if coeff > 0 else -1
            else:
                converted_efm[r_id] = coeff
        if math.log(i) == int_size:
            cur_efm_part = next(binary_efm_iterable)
            i = 1
        else:
            i <<= 1
        return i, cur_efm_part

    i = 1
    for r_id in r_ids:
        i, cur_efm_part = process(r_id, cur_efm_part, i)
        if r_id in reversible_r_ids:
            i, cur_efm_part = process(r_id, cur_efm_part, i, True)
    return converted_efm


def efm2binary(efm, r_ids, reversible_r_ids, int_size=None):
    """
    Converts an EFM representation {r_id: coefficient} into a tuple (binary_representation, non_zero_coefficients).

    A binary representation of an EFM is a list of integers whose binary representations
    correspond to the reactions that are active in the EFM: if the reaction is active,
    the corresponding bit is set to 1.
    If the total number of reactions in the model is larger that the number of bits in an int, several ints are used.

    Example: For a model containing reactions r1, r2(reversible), r3(reversible), r4, r5,
    a EFM: 3 r1, -2 r2, 1 r3, 1 r5 would be represented as [77], as the binary representation of 77 is '1001101'
    that corresponds to '1 r5, 0 r4, 0 -r3, 1 r3, 1 -r2, 0 r2, 1 r1', where -ri is the reversed version of ri.
    The non_zero_coefficients for this EFM are [3, -2, 1, 1].

    :param efm: a EFM represented as a dictionary {r_id: coefficient}.

    :param r_ids: ordered collection of reaction ids (strings).

    :param reversible_r_ids: set of ids of reversible reactions (strings).

    :param int_size: the maximal number of bits in an int, can be calculated as math.log(sys.maxint) / math.log(2).

    :return: a EFM represented as a tuple: binary_representation, non_zero_coefficients.
    """
    if not int_size:
        int_size = get_int_size()
    bit_efm = []
    bit_efm_cur = 0
    i = 1
    coefficients = []

    def advance(i, bit_efm_cur):
        if math.log(i) == int_size:
            bit_efm.append(bit_efm_cur)
            bit_efm_cur = 0
            i = 1
        else:
            i <<= 1
        return bit_efm_cur, i

    for r_id in r_ids:
        if r_id in efm:
            coeff = efm[r_id]
        else:
            coeff = 0
        if coeff > 0:
            bit_efm_cur |= i
            coefficients.append(coeff)
        bit_efm_cur, i = advance(i, bit_efm_cur)
        if r_id in reversible_r_ids:
            if coeff < 0:
                bit_efm_cur |= i
                coefficients.append(coeff)
            bit_efm_cur, i = advance(i, bit_efm_cur)
    bit_efm.append(bit_efm_cur)
    return tuple(bit_efm), coefficients


def get_int_size():
    """
    Calculates the maximal number of bits in an int:
    math.log(sys.maxint) / math.log(2).

    :return: int, the maximal number of bits in an int.
    """
    return math.log(sys.maxint) / math.log(2)


def get_binary_efm_length(binary_efm):
    """
    Returns the length of a EFM represented in a binary form (number of active reactions).

    A binary representation of an EFM is a list of integers whose binary representations
    correspond to the reactions that are active in the EFM: if the reaction is active,
    the corresponding bit is set to 1.
    If the total number of reactions in the model is larger that the number of bits in an int, several ints are used.

    Example: For a model containing reactions r1, r2(reversible), r3(reversible), r4, r5,
    a EFM: 3 r1, -2 r2, 1 r3, 1 r5 would be represented as [77], as the binary representation of 77 is '1001101'
    that corresponds to '1 r5, 0 r4, 0 -r3, 1 r3, 1 -r2, 0 r2, 1 r1', where -ri is the reversed version of ri.

    :param binary_efm: list of integers -- an EFM represented in a binary form.
    :return: int, number of active reactions in the EFM.
    """
    return sum(gmpy.popcount(it) for it in binary_efm)


def compute_efms(sbml, directory, em_number, r_id, rev, tree_efm_path, r_id2rev_2threshold=None, threshold=0.0,
                 r_ids=None):
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
    :return:efms: list of EFMs in binary representation: (binary_efm, non_zero_coefficients);
            r_ids: concerned reaction ids;


    A binary representation of an EFM is a list of integers whose binary representations
    correspond to the reactions that are active in the EFM: if the reaction is active,
    the corresponding bit is set to 1.
    If the total number of reactions in the model is larger that the number of bits in an int, several ints are used.

    Example: For a model containing reactions r1, r2(reversible), r3(reversible), r4, r5,
    a EFM: 3 r1, -2 r2, 1 r3, 1 r5 would be represented as [77], as the binary representation of 77 is '1001101'
    that corresponds to '1 r5, 0 r4, 0 -r3, 1 r3, 1 -r2, 0 r2, 1 r1', where -ri is the reversed version of ri.
    The non_zero_coefficients for this EFM are [3, -2, 1, 1].

    :raise ValueError: if the reaction of interest was not found in the model.
    """
    if not os.path.exists(tree_efm_path):
        raise ValueError("TreeEFM runner is not found at %s" % tree_efm_path)
    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()
    # Compute stoichiometric matrix
    st_matrix_file = "%s/st_matrix.txt" % directory
    s_id2i, r_id2i, rev_r_id2i = stoichiometric_matrix(model, st_matrix_file)
    logging.info("stoichiometric matrix saved to %s" % st_matrix_file)
    # Figure out in which reaction we are interested in
    if (rev and r_id not in rev_r_id2i) or (not rev and r_id not in r_id2i):
        raise ValueError("R%seaction with id %s is not found in the model" % ('eversed r' if rev else '', r_id))
    i = rev_r_id2i[r_id] if rev else r_id2i[r_id]
    # Compute EFMs using TreeEFM software (Pey et al. 2014, PMID: 25380956)
    em_file = "%s/FV-EM.dat" % directory
    os.system("%s -r %d -i %s -l EM -e %d -o %s" % (tree_efm_path, i, st_matrix_file, em_number, directory))
    os.system("%s -b %s" % (tree_efm_path, em_file))
    em_file = "%s.txt" % em_file
    logging.info("elementary modes saved to %s" % em_file)
    # Filter EFMs so that only those that don't include the reaction in opposite directions are left.
    # If r_id2rev are specified, filter EFMs to leave only those that include these reactions in these directions.
    em_file_filtered = "%s/FV-EM_filtered.dat.txt" % directory
    efms = filter_efms(em_file, r_id2i, rev_r_id2i, em_file_filtered, r_id2rev_2threshold, zero_threshold=threshold,
                      r_ids=r_ids)
    logging.info(
        "%d elementary modes corresponding to reactions of interest saved to %s" % (len(efms), em_file_filtered))
    efms = sorted(efms, key=lambda (binary_efm, coeffs): len(coeffs))
    return efms, r_ids if r_ids else sorted(set(r_id2i.iterkeys()) | set(rev_r_id2i.iterkeys()))


def filter_efms(in_path, r_id2i, rev_r_id2i, out_path, r_id2rev_2threshold, zero_threshold=0.0, r_ids=None):
    i2r_id = {i: (r_id, False) for (r_id, i) in r_id2i.iteritems()}
    i2r_id.update({i: (r_id, True) for (r_id, i) in rev_r_id2i.iteritems()})
    all_r_ids = r_ids if r_ids else sorted(set(r_id2i.iterkeys()) | set(rev_r_id2i.iterkeys()))
    rev_r_ids = set(rev_r_id2i.iterkeys())
    int_size = math.log(sys.maxint) / math.log(2)
    processed = set()
    efms = []
    round_value = lambda v: round(float(v), PRECISION)
    with open(out_path, 'w+') as out_f:
        with open(in_path, 'r') as in_f:
            for line in in_f:
                values = line.replace("\n", "").strip().split(" ")
                efm = {}
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
                    if r_id in efm:
                        bad_em = True
                        break
                    efm[r_id] = v

                # The same reaction participates in different directions
                # => don't want such an EFM
                if bad_em:
                    continue

                # Check if all the reactions that are required are there
                bad_em = False
                if r_id2rev_2threshold:
                    for (r_id2rev, present_reaction_threshold) in r_id2rev_2threshold:
                        present_r_ids = set()
                        for (r_id, rev) in r_id2rev.iteritems():
                            if r_id in efm and (rev is None or rev == (efm[r_id] < 0)):
                                present_r_ids.add(r_id)
                        if len(present_r_ids) * 1.0 / len(r_id2rev) < present_reaction_threshold:
                            bad_em = True
                            break
                    if bad_em:
                        continue

                out_f.write(line)

                if r_ids:
                    efm = {r_id: val for (r_id, val) in efm.iteritems() if r_id in r_ids}
                    # em_support = tuple(sorted((r_id, val > 0) for (r_id, val) in em.iteritems()))
                em_support, coefficients = efm2binary(efm, all_r_ids, rev_r_ids, int_size)
                if not r_ids or em_support not in processed:
                    efms.append((em_support, coefficients))
                    if r_ids:
                        processed.add(em_support)
    return efms


def classify_efms_with_acom(efms, min_motif_length, r_ids, neighbour_threshold=None):
    """
    Classifies EFMs to find common motifs
    (using ACoM method [Peres et al. 2011, doi:10.1016/j.biosystems.2010.12.001]).

    :param efms: dictionary {r_id: coefficient] of EFMs.
    :param min_motif_length:int, minimal motif length to be considered for ACoM classification (see Peres et al. 2011)
    :param r_ids: collection of reaction ids.
    :param neighbour_threshold: int, at least how many common elements two EFMs should have
    to be considered as neighbours in ACoM classification (see Peres et al. 2011).
    If not specified, then is set to the mean of all the values of the resemblance matrix.
    :return: list of clusters: [(motif, list of binary EFMs that contain this motif)]
    and a list of outliers: binary EFMs that were not clustered.
    """
    r_ids = sorted(r_ids)
    # Convert EFMs to binary EFMs
    binary_ems = ems2binary(efms, r_ids)
    logging.info(
        "elementary modes converted to %d binary vectors" % len(binary_ems))
    # Classify binary EFMs using ACoM method (Peres et al. 2011, doi:10.1016/j.biosystems.2010.12.001)
    clusters, outliers = classify(binary_ems, min_motif_length, neighbour_threshold)
    logging.info("---------CLUSTERS (%d)-------" % len(clusters))
    clusters = {tuple(sorted((value, r_ids[index]) for (value, index) in motif)): elements for (motif, elements) in
                clusters.iteritems()}
    for motif, elements in clusters.iteritems():
        logging.info(
            "%d %d %s %s" % (len(motif), len(elements), motif, elements))
    logging.info("---------OUTLIERS (%d)-------" % len(outliers))
    logging.info(sorted(outliers))

    return clusters, outliers
