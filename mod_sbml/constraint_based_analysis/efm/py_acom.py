from collections import defaultdict
import logging
from itertools import islice

__author__ = 'anna'


def get_motif(e, l=None):
    if not l:
        l = len(e)
    return {(value, index) for (value, index) in zip(e, xrange(l)) if value != 0}


def calculate_resemblance_matrix(ems, neighbour_threshold):
    resemblance_matrix = []
    s, l = 0, 0
    i = 0
    e_len = len(ems[0])
    i2motif = {}
    for em in ems:
        motif = get_motif(em, e_len)
        i2motif[i] = motif
        values = [len(motif & i2motif[j]) for j in xrange(0, i + 1)]
        s += sum(values)
        l += len(values)
        resemblance_matrix.append(values)
        i += 1
    if not neighbour_threshold:
        neighbour_threshold = round(s * 1.0 / l)
    logging.info("using %d as neighbour threshold" % neighbour_threshold)
    return i2motif, resemblance_matrix, neighbour_threshold


def calculate_clusters(ems, i2motif, resemblance_matrix, min_motif_length, neighbour_threshold):
    clusters, outliers = defaultdict(set), set()
    ems_len = len(ems)
    outliers = set(xrange(ems_len))
    for i in xrange(ems_len):
        stop_it = False
        neighbours = {i}
        motif = set(i2motif[i])
        j = 0
        for j_resemblance in resemblance_matrix[i]:
            if j_resemblance < neighbour_threshold:
                continue
            motif &= i2motif[j]
            if len(motif) < min_motif_length:
                stop_it = True
                break
            neighbours.add(j)
            j += 1
        if stop_it:
            outliers.add(i)
            continue
        j = i + 1
        for row in islice(resemblance_matrix, i + 1, None):
            if row[i] < neighbour_threshold:
                continue
            motif &= i2motif[j]
            if len(motif) < min_motif_length:
                stop_it = True
                break
            neighbours.add(j)
            j += 1
        if stop_it:
            outliers.add(i)
            continue
        if len(neighbours) > 1:
            clusters[tuple(sorted(motif))] |= neighbours
            outliers -= neighbours
    logging.info("found %d clusters" % len(clusters))
    logging.info("going to cluster %d outliers" % len(outliers))
    return clusters, outliers


def merge_similar_clusters(clusters, min_motif_length):
    changed = True
    while changed:
        changed = False
        old_keys = clusters.keys()
        i = 0
        for motif in old_keys:
            other_keys = list(islice(old_keys, i + 1, None))
            for motif2 in other_keys:
                if motif2 == motif:
                    continue
                intersection = set(motif) & set(motif2)
                if len(intersection) >= min_motif_length:
                    intersection = tuple(sorted(intersection))
                    if intersection in clusters:
                        clusters[intersection] |= clusters[motif] | clusters[motif2]
                    else:
                        clusters[intersection] = clusters[motif] | clusters[motif2]
                    if intersection != motif:
                        del clusters[motif]
                    if intersection != motif2:
                        del clusters[motif2]
                    changed = True
                    motif = intersection
            i += 1
            if changed:
                break


def classify_outliers(clusters, outliers, i2motif, min_motif_length):
    while clusters:
        clustered_ems = set()
        for i in outliers:
            old_keys = clusters.keys()
            for motif in old_keys:
                c_motif = set(motif) & i2motif[i]
                if len(c_motif) >= min_motif_length:
                    if len(motif) != len(c_motif):
                        c_motif = tuple(sorted(c_motif))
                        elements = clusters[motif]
                        del clusters[motif]
                        clusters[c_motif] |= elements | {i}
                    else:
                        clusters[motif].add(i)
                    clustered_ems.add(i)
        outliers -= clustered_ems
        if not clustered_ems or not outliers:
            break
    return outliers


def classify(ems, min_motif_length, neighbour_threshold=None):
    """
    Classifies binary elementary flux modes (EFMs) using ACoM method:
    Peres et al. 2011 (doi:10.1016/j.biosystems.2010.12.001).
    :param ems: a list of binary EFMa,
    i.e. for each binary EFM, EFM[i] = 1 | -1 | 0 if coefficient(reaction i) > 0 | < 0 | == 0.
    :param min_motif_length: int, minimal motif length to be considered.
    :param neighbour_threshold: int, at least how many common elements two EFMs should have
    to be considered as neighbours (see Peres et al. 2011).
    If not specified, then is set to the mean of all the values of the resemblance matrix.
    :return: list of clusters: [(motif, list of binary EFMs that contain this motif)]
    and a list of outliers: binary EFMs that were not clustered.
    """
    if not ems:
        return None, None
    i2motif, resemblance_matrix, neighbour_threshold = calculate_resemblance_matrix(ems, neighbour_threshold)
    clusters, outliers = calculate_clusters(ems, i2motif, resemblance_matrix, min_motif_length, neighbour_threshold)
    merge_similar_clusters(clusters, min_motif_length)
    outliers = classify_outliers(clusters, outliers, i2motif, min_motif_length)
    return clusters, outliers


def classify_efms_with_acom(efms, r_ids, rev_r_ids, min_motif_length, neighbour_threshold=None):
    """
    Classifies EFMs to find common motifs
    (using ACoM method [Peres et al. 2011, doi:10.1016/j.biosystems.2010.12.001]).

    :param efms: a list of EFMs in binary form.

    A binary representation of an EFM is a list of integers whose binary representations
    correspond to the reactions that are active in the EFM: if the reaction is active,
    the corresponding bit is set to 1.
    If the total number of reactions in the model is larger that the number of bits in an int, several ints are used.

    Example: For a model containing reactions r1, r2(reversible), r3(reversible), r4, r5,
    a EFM: 3 r1, -2 r2, 1 r3, 1 r5 would be represented as [77], as the binary representation of 77 is '1001101'
    that corresponds to '1 r5, 0 r4, 0 -r3, 1 r3, 1 -r2, 0 r2, 1 r1', where -ri is the reversed version of ri.

    :param r_ids: ordered collection of reaction ids (strings).

    :param rev_r_ids: set of ids of reversible reactions (strings).

    :param min_motif_length:int, minimal motif length to be considered for ACoM classification (see Peres et al. 2011)

    :param neighbour_threshold: int, at least how many common elements two EFMs should have
    to be considered as neighbours in ACoM classification (see Peres et al. 2011).
    If not specified, then is set to the mean of all the values of the resemblance matrix.

    :return: list of clusters: [(motif, list of binary EFMs that contain this motif)]
    and a list of outliers: binary EFMs that were not clustered.
    """
    r_ids = sorted(r_ids)
    # Convert EFMs to binary EFMs
    binary_ems = [efm.binary_efm for efm in efms]
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

