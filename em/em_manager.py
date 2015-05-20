import logging
import os
import datetime

import libsbml
import openpyxl

from em.acom import classify
from em.stoichiometry_manager import stoichiometric_matrix, ems2binary
from gibbs.reaction_boundary_manager import get_bounds
from sbml.sbml_manager import get_kegg_r_id, get_gene_association, get_r_comps, submodel
from serialization.serialization_manager import get_sbml_r_formula
from serialization.xlsx_helper import save_data, BASIC_STYLE


basic_r_style = lambda r_id: BASIC_STYLE

__author__ = 'anna'


ZERO_THRESHOLD = 1e-9


def compute_efms(sbml, directory, em_number, r_id, rev, tree_efm_path, r_id2rev_2threshold=None, output_file=None,
                 over_expressed_r_ids=set(), under_expressed_r_ids=set(), threshold=0.0, r_ids=None,
                 r_id2style=basic_r_style):
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
    :param output_file: string, if specified, the EFMs will be saved to this xlsx file:
    each EFM in a separate sheet called EFM_<EFM_number>_<number_of_participating_reactions>,
    the sheet contains the information about the reactions in the EFM.
    :return:efms: dictionary r_id: coefficient (only contains values with non-zero coefficients);
            r_ids: concerned reaction ids;
    :raise ValueError: if the reaction of interest was not found in the model.
    """
    if not os.path.exists(tree_efm_path):
        raise ValueError("TreeEFM runner is not found at %s" % tree_efm_path)
    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()
    # Compute stoichiometric matrix
    st_matrix_file = "%s/st_matrix.txt" % directory
    s_id2i, r_id2i = stoichiometric_matrix(model, st_matrix_file)
    logging.info("%s\nstoichiometric matrix saved to %s" % (datetime.datetime.now().time(), st_matrix_file))
    # Figure out in which reaction we are interested in
    if r_id not in r_id2i:
        raise ValueError("Reaction with id %s is not found in the model" % r_id)
    i = r_id2i[r_id][1] if rev else r_id2i[r_id][0]
    # Compute EFMs using TreeEFM software (Pey et al. 2014, PMID: 25380956)
    em_file = "%s/FV-EM.dat" % directory
    os.system("%s -r %d -i %s -l EM -e %d -o %s" % (tree_efm_path, i, st_matrix_file, em_number, directory))
    os.system("%s -b %s" % (tree_efm_path, em_file))
    em_file = "%s.txt" % em_file
    logging.info("%s\nelementary modes saved to %s" % (datetime.datetime.now().time(), em_file))
    # Filter EFMs so that only those that don't include the reaction in opposite directions are left.
    # If r_id2rev are specified, filter EFMs to leave only those that include these reactions in these directions.
    em_file_filtered = "%s/FV-EM_filtered.dat.txt" % directory
    efms = filter_ems(em_file, r_id2i, em_file_filtered, r_id2rev_2threshold, threshold=threshold, r_ids=r_ids)
    logging.info("%s\n%d elementary modes corresponding to reactions of interest saved to %s" % (
        datetime.datetime.now().time(), len(efms), em_file_filtered))
    # em_file = em_file_filtered
    # efms = sorted(format_ems(em_file, r_id2i, threshold, r_ids_to_keep=r_id_of_interest), key=len)
    efms = sorted(efms, key=len)
    # Save the result to file
    if output_file:
        serialize_efms(sbml, efms, output_file, over_expressed_r_ids=over_expressed_r_ids,
                       under_expressed_r_ids=under_expressed_r_ids, r_id2style=r_id2style)

    return efms, r_ids if r_ids else r_id2i.keys()


def filter_ems(in_path, r_id2i, out_path, r_id2rev_2threshold, threshold=0.0, r_ids=None):
    i2r_id = {}
    for r_id, (i, j) in r_id2i.iteritems():
        if i:
            i2r_id[i] = (r_id, False)
        if j:
            i2r_id[j] = (r_id, True)
    num = 0
    processed = set()
    ems = []
    with open(out_path, 'w+') as out_f:
        with open(in_path, 'r') as in_f:
            for line in in_f:
                values = line.replace("\n", "").strip().split(" ")
                em = {i2r_id[i] for (v, i) in zip(values, xrange(1, len(values) + 1)) if round(float(v), 6)}
                bad_em = False
                for (r_id, rev) in em:
                    if (r_id, not rev) in em:
                        bad_em = True
                        break
                if bad_em:
                    continue
                if r_id2rev_2threshold:
                    for (r_id2rev, present_reaction_threshold) in r_id2rev_2threshold:
                        present_r_ids = set()
                        em_r_ids = {it[0] for it in em}
                        for (r_id, rev) in r_id2rev.iteritems():
                            if (r_id, rev) in em or (rev is None and r_id in em_r_ids):
                                present_r_ids.add(r_id)
                        if len(present_r_ids) * 1.0 / len(r_id2rev) < present_reaction_threshold:
                            bad_em = True
                            break
                    if bad_em:
                        continue
                out_f.write(line)
                em = {i2r_id[i][0]: (-round(float(v), 6) if i2r_id[i][1] else round(float(v), 6))
                      for (v, i) in zip(values, xrange(1, len(values) + 1))
                      if abs(round(float(v), 6)) > threshold}
                if r_ids:
                    em = {r_id: val for (r_id, val) in em.items() if r_id in r_ids}
                    em_support = tuple(sorted((r_id, val > 0) for (r_id, val) in em.iteritems()))
                if not r_ids or em_support not in processed:
                    ems.append(em)
                    if r_ids:
                        processed.add(em_support)
                    num += 1
    return ems


def classify_efms(efms, min_motif_length, r_ids, neighbour_threshold=None, output_file=None, sbml=None,
                  r_id2style=basic_r_style):
    """
    Classifies EFMs to find common motifs
    (using ACoM method [Peres et al. 2011, doi:10.1016/j.biosystems.2010.12.001]).

    :param efms: dictionary {r_id: coefficient] of EFMs.
    :param min_motif_length:int, minimal motif length to be considered for ACoM classification (see Peres et al. 2011)
    :param r_ids: collection of reaction ids.
    :param neighbour_threshold: int, at least how many common elements two EFMs should have
    to be considered as neighbours in ACoM classification (see Peres et al. 2011).
    If not specified, then is set to the mean of all the values of the resemblance matrix.
    :param output_file: string, if specified, the classes of EFMs will be saved to this xlsx file:
    each motif in a separate sheet called Motif_<motif_number>_<motif_length>_<number_of_EFMs_that_contain_this_motif>,
    the sheet contains the information about the reactions in the motif.
    :param sbml: string, path to the SBML file with the model
    (needed for saving EFMs to output file, otherwise can be None).
    :return: list of clusters: [(motif, list of binary EFMs that contain this motif)]
    and a list of outliers: binary EFMs that were not clustered.
    """
    r_ids = sorted(r_ids)
    # Convert EFMs to binary EFMs
    binary_ems = ems2binary(efms, r_ids)
    logging.info(
        "%s\nelementary modes converted to %d binary vectors" % (datetime.datetime.now().time(), len(binary_ems)))
    # Classify binary EFMs using ACoM method (Peres et al. 2011, doi:10.1016/j.biosystems.2010.12.001)
    clusters, outliers = classify(binary_ems, min_motif_length, neighbour_threshold)
    logging.info('''%s
    classification done
    ---------CLUSTERS (%d)-------
    ''' % (datetime.datetime.now().time(), len(clusters)))
    clusters = {tuple(sorted((value, r_ids[index]) for (value, index) in motif)): elements for (motif, elements) in
                clusters.iteritems()}
    for motif, elements in clusters.iteritems():
        logging.info(
            "%d %d %s %s" % (len(motif), len(elements), motif, elements))
    logging.info("---------OUTLIERS (%d)-------" % len(outliers))
    logging.info(outliers)

    # Save the result to file
    if output_file and sbml:
        doc = libsbml.SBMLReader().readSBML(sbml)
        model = doc.getModel()
        wb = openpyxl.Workbook()
        i = 0
        for motif in clusters.iterkeys():
            save_data(["Id", "Name", "Formula", "Lower bound", "Upper bound", "Direction"],
                      [[r.id, r.name, get_sbml_r_formula(model, r, False),
                        get_bounds(r)[0], get_bounds(r)[1], direction]
                       for (direction, r) in ((direction, model.getReaction(r_id)) for (direction, r_id) in motif)],
                      wb, "Motif_%d_%d_%d.xlsx" % (i, len(motif), len(clusters[motif])), i,
                      styles=[r_id2style(r_id) for (_, r_id) in motif])
            i += 1
        wb.save(output_file)

    return clusters, outliers


def analyse_ems(sbml, directory, tree_efm_path, r_id, min_motif_length, neighbour_threshold=None, rev=False,
                em_number=100, output_efm_file=None, output_motif_file=None, r_id2rev_2threshold=None,
                over_expressed_r_ids=set(), under_expressed_r_ids=set(), threshold=0.0, r_ids=None,
                r_id2style=basic_r_style):
    """
    Analyses elementary flux modes (EFMs) in a given SBML (see http://sbml.org) model:
    1. computes the EFMs that contain a reaction of interest
    (using TreeEFM software [Pey et al. 2014, PMID: 25380956]);
    2. classifies them to find common motifs
    (using ACoM method [Peres et al. 2011, doi:10.1016/j.biosystems.2010.12.001]).

    :param sbml: string, path to the SBML file with the model.
    :param directory: string, directory where to store the results, such as stoichiometric matrix, EFMs.
    :param tree_efm_path: string,path to the executable of TreeEFM software [Pey et al. 2014, PMID: 25380956].
    :param r_id:string, id of the reaction of interest.
    :param min_motif_length:int, minimal motif length to be considered for ACoM classification (see Peres et al. 2011)
    :param neighbour_threshold: int, at least how many common elements two EFMs should have
    to be considered as neighbours in ACoM classification (see Peres et al. 2011).
    If not specified, then is set to the mean of all the values of the resemblance matrix.
    :param rev: boolean, if the reaction of interest should be considered in the opposite direction.
    :param em_number: int, number of EFMs to compute with TreeEFM software [Pey et al. 2014, PMID: 25380956].
    :param output_efm_file: string, if specified, the EFMs will be saved to this xlsx file:
    each EFM in a separate sheet called EFM_<EFM_number>_<number_of_participating_reactions>,
    the sheet contains the information about the reactions in the EFM.
    :param output_motif_file: string, if specified, the classes of EFMs will be saved to this xlsx file:
    each motif in a separate sheet called Motif_<motif_number>_<motif_length>_<number_of_EFMs_that_contain_this_motif>,
    the sheet contains the information about the reactions in the motif.
    :param r_id2rev_2threshold: set of strings in the form {(r_id_0, reversed_0), (r_id_1, reversed_1), ...},
    if specified, only EFMs that contain all (if reaction_op == (REACTION_OPERATION_AND))
    or any (if reaction_op == (REACTION_OPERATION_OR) of the reaction ids in specified directions
    from this set will be returned.
    :return: list of clusters: [(motif, list of binary EFMs that contain this motif)];
    list of outliers: binary EFMs that were not clustered;
    dictionary r_id: coefficient (only contains values with non-zero coefficients).
    :raise ValueError: if the reaction of interest was not found in the model.
    """
    efms, r_ids = compute_efms(sbml, directory, em_number, r_id, rev, tree_efm_path, r_id2rev_2threshold,
                               output_efm_file, over_expressed_r_ids=over_expressed_r_ids,
                               under_expressed_r_ids=under_expressed_r_ids, threshold=threshold, r_ids=r_ids,
                               r_id2style=r_id2style)

    clusters, outliers = classify_efms(efms, min_motif_length, r_ids, neighbour_threshold, output_motif_file, sbml,
                                       r_id2style=r_id2style)

    return clusters, outliers, efms


def serialize_efms(sbml, efms, path, r_id2style=basic_r_style,
                   over_expressed_r_ids=set(), under_expressed_r_ids=set()):
    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()
    wb = openpyxl.Workbook()
    i = 0
    for r_id2coefficients in sorted(efms, key=lambda r_id2coeffs: (
            -len(set(r_id2coeffs.iterkeys()) & over_expressed_r_ids),
            len(set(r_id2coeffs.iterkeys()) & under_expressed_r_ids),
            len(r_id2coeffs))):
        data, styles = [], []
        for r_id in sorted(r_id2coefficients.iterkeys(), key=lambda r_id: (sorted(get_r_comps(r_id, model)), r_id)):
            r = model.getReaction(r_id)
            if not r:
                raise ValueError('Reaction with id %s was not found in the model %s' % (r_id, model.getId()))
            lb, ub = get_bounds(r)
            comps = ", ".join(sorted((model.getCompartment(c_id).name for c_id in get_r_comps(r_id, model))))
            data.append([r.id, r.name, get_sbml_r_formula(model, r, False), get_kegg_r_id(r), get_gene_association(r),
                         lb, ub, ','.join(sorted(get_r_comps(r.id, model))), r_id2coefficients[r_id], comps])
            styles.append(r_id2style(r.id))
        r_ids = set(r_id2coefficients.iterkeys())
        save_data(["Id", "Name", "Formula", "Kegg", "Genes", "Low. B.", "Upp. B.", "Compartments", "Coefficients"],
                  data=data, ws_name="EM_%d_(%d)_%d_%d" % (i + 1, len(r_id2coefficients),
                                                           len(r_ids & over_expressed_r_ids),
                                                           len(r_ids & under_expressed_r_ids)),
                  wb=wb, ws_index=i, styles=styles)
        i += 1
    wb.save(path)


def perform_efma(sbml, in_r_id, in_r_reversed, out_r_id2rev_2threshold, em_number,
                   min_motif_len, model_dir, neighbour_threshold=None, output_motif_file=None, output_efm_file=None,
                   motif_sbml_prefix=None, efm_sbml_prefix=None,
                   over_expressed_r_ids=set(), under_expressed_r_ids=set(), threshold=ZERO_THRESHOLD, r_ids=None,
                   r_id2style=basic_r_style, tree_efm_path="/home/anna/Applications/TreeEFM/tool/TreeEFMseq"):
    clusters, outliers, efms = analyse_ems(sbml, model_dir, em_number=em_number,
                                           r_id=in_r_id, rev=in_r_reversed,
                                           tree_efm_path=tree_efm_path,
                                           min_motif_length=min_motif_len, neighbour_threshold=neighbour_threshold,
                                           r_id2rev_2threshold=out_r_id2rev_2threshold,
                                           output_motif_file=output_motif_file, output_efm_file=output_efm_file,
                                           over_expressed_r_ids=over_expressed_r_ids,
                                           under_expressed_r_ids=under_expressed_r_ids,
                                           threshold=threshold, r_ids=r_ids, r_id2style=r_id2style)

    if motif_sbml_prefix is not None:
        i = 1
        for motif in clusters.iterkeys():
            doc = libsbml.SBMLReader().readSBML(sbml)
            model = doc.getModel()
            submodel({r_id for (direction, r_id) in motif}, model)
            model.setId('%s_Motif_%d' % (model.getId(), i))
            model.setName('%s_Motif_%d' % (model.getName(), i))
            for (direction, r_id) in motif:
                r = model.getReaction(r_id)
                r.setName('%s%s' % ('-' if -1 == direction else '', r.getId()))
            libsbml.SBMLWriter().writeSBMLToFile(doc, '%s_%d_%d_%d.xml'
                                                 % (motif_sbml_prefix, i, len(clusters[motif]), len(motif)))
            i += 1

    if efm_sbml_prefix is not None:
        i = 1
        for r_id2coeff in efms:
            doc = libsbml.SBMLReader().readSBML(sbml)
            model = doc.getModel()
            submodel(set(r_id2coeff.iterkeys()), model)
            model.setId('%s_EM_%d_%d' % (model.getId(), i, len(r_id2coeff)))
            model.setName('%s_EM_%d_%d' % (model.getName(), i, len(r_id2coeff)))
            for (r_id, coeff) in r_id2coeff.iteritems():
                r = model.getReaction(r_id)
                r.setName('%s %g' % (r.getName(), coeff))
            libsbml.SBMLWriter().writeSBMLToFile(doc, '%s_%d_%d.xml' % (efm_sbml_prefix, i, len(r_id2coeff)))
            i += 1
