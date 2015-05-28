from collections import defaultdict
import logging
import math
import sys

import libsbml
import openpyxl

from em.em_manager import binary2efm, get_binary_efm_length, get_int_size
from sbml.sbml_manager import get_r_comps, submodel, get_kegg_r_id
from serialization.serialization_manager import get_sbml_r_formula
from serialization.xlsx_helper import BASIC_STYLE, save_data

SIMPLE_PATTERN_SORTER = lambda p_id: -p_id

__author__ = 'anna'

basic_r_style = lambda r_id: BASIC_STYLE

SORT_BY_PATTERN_LENGTH = 1
SORT_BY_EFM_NUMBER = 0
SORT_BY_WEIGHTED_PRODUCT_OF_EFM_NUMBER_AND_PATTERN_LENGTH = 2


def get_pattern_sorter(id2pattern, p_id2efm_ids, min_pattern_len, min_efm_num, sort=SORT_BY_EFM_NUMBER):
    if not id2pattern:
        return SIMPLE_PATTERN_SORTER
    if SORT_BY_PATTERN_LENGTH == sort:
        return lambda p_id: -get_binary_efm_length(id2pattern[p_id])
    elif SORT_BY_EFM_NUMBER == sort:
        return lambda p_id: -len(p_id2efm_ids[p_id])
    elif SORT_BY_WEIGHTED_PRODUCT_OF_EFM_NUMBER_AND_PATTERN_LENGTH == sort:
        max_pattern_len = max(get_binary_efm_length(pattern) for pattern in id2pattern.itervalues())
        max_efm_number = max(len(efm_ids) for efm_ids in p_id2efm_ids.itervalues())
        return lambda p_id: \
            -(1.0 * (get_binary_efm_length(id2pattern[p_id]) - min_pattern_len) / (max_pattern_len - min_pattern_len)
              if max_pattern_len != min_pattern_len else 1) * \
            ((1.0 * len(p_id2efm_ids[p_id]) - min_efm_num) / (max_efm_number - min_efm_num)
             if max_efm_number != min_efm_num else 1)
    else:
        logging.error('Unknown pattern sorting function %s. Will sort by pattern id.' % sort)
        return SIMPLE_PATTERN_SORTER


def serialize_efms_xslx(sbml, efms, r_ids, reversible_r_ids, path, r_id2style=basic_r_style,
                        over_expressed_r_ids=set(), under_expressed_r_ids=set()):
    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()
    wb = openpyxl.Workbook()
    int_size = math.log(sys.maxint) / math.log(2)
    i = 0
    for (binary_efm, coefficients) in sorted(efms, key=lambda (_, coefficients): len(coefficients)):
        r_id2coefficient = binary2efm(binary_efm, r_ids, reversible_r_ids, int_size, coefficients)
        data, styles = [], []
        for r_id in sorted(r_id2coefficient.iterkeys(), key=lambda r_id: (sorted(get_r_comps(r_id, model)), r_id)):
            r = model.getReaction(r_id)
            if not r:
                raise ValueError('Reaction with id %s was not found in the model %s' % (r_id, model.getId()))
            comps = ", ".join(sorted((c.name if c.name else c.id for c in
                                      (model.getCompartment(c_id) for c_id in get_r_comps(r_id, model)))))
            data.append([r_id2coefficient[r_id], r.id, r.name, get_sbml_r_formula(model, r, False), get_kegg_r_id(r),
                         comps])
            styles.append(r_id2style(r.id))
        efm_r_ids = set(r_id2coefficient.iterkeys())
        save_data(["Coefficient", "Id", "Name", "Formula", "Kegg", "Compartments"],
                  data=data, ws_name="EM_%d_(%d)_%d_%d" % (i + 1, len(r_id2coefficient),
                                                           len(efm_r_ids & over_expressed_r_ids),
                                                           len(efm_r_ids & under_expressed_r_ids)),
                  wb=wb, ws_index=i, styles=styles)
        i += 1
    wb.save(path)


def serialize_efms_txt(sbml, efms, r_ids, reversible_r_ids, path):
    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()
    int_size = math.log(sys.maxint) / math.log(2)
    with open(path, 'w+') as f:
        for (binary_efm, coefficients) in sorted(efms, key=lambda (_, coefficients): len(coefficients)):
            r_id2coefficient = binary2efm(binary_efm, r_ids, reversible_r_ids, int_size, coefficients)
            f.write(', '.join('%g %s' % (r_id2coefficient[r_id], r_id)
                              for r_id in sorted(r_id2coefficient.iterkeys(),
                                                 key=lambda r_id: (sorted(get_r_comps(r_id, model)), r_id))) + '\n')


def efm2sbml(id2efm, get_name_suffix, name_prefix, directory, sbml, r_ids, rev_r_ids, binary):
    int_size = get_int_size()
    for efm_id, efm in id2efm.iteritems():
        if isinstance(efm, tuple) and len(efm) == 2 and isinstance(efm[1], list):
            efm, coeffs = efm
        else:
            coeffs = None
        r_id2coeff = binary2efm(efm, r_ids, rev_r_ids, int_size, coeffs, binary)
        doc = libsbml.SBMLReader().readSBML(sbml)
        model = doc.getModel()
        submodel(set(r_id2coeff.iterkeys()), model)
        suffix = get_name_suffix(efm_id)
        model.setId('%s_%s_%d_%s' % (model.getId(), name_prefix, efm_id, suffix))
        model.setName('%s_%s_%d_%s' % (model.getName(), name_prefix, efm_id, suffix))
        for (r_id, coeff) in r_id2coeff.iteritems():
            r = model.getReaction(r_id)
            r_name = r.getName()
            r.setName('%g %s' % (coeff, r_name) if abs(coeff) != 1 else (r_name if coeff > 0 else ('-' + r_name)))
        libsbml.SBMLWriter().writeSBMLToFile(doc, '%s/%s_%d_%s.xml' % (directory, name_prefix, efm_id, suffix))


def serialize_acom_patterns(clusters, output_file, sbml, r_id2style=basic_r_style):
    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()
    wb = openpyxl.Workbook()
    i = 0
    for motif in clusters.iterkeys():
        data, styles = [], []
        for (direction, r_id) in motif:
            r = model.getReaction(r_id)
            if not r:
                raise ValueError('Reaction with id %s was not found in the model %s' % (r_id, model.getId()))
            data.append([r.id, r.name, get_sbml_r_formula(model, r, False), direction])
            styles.append(r_id2style(r.id))
        save_data(["Id", "Name", "Formula", "Direction"],
                  data, wb, "Motif_%d_%d_%d.xlsx" % (i, len(motif), len(clusters[motif])), i,
                  styles=styles)
        i += 1
    wb.save(output_file)


def serialize_patterns(p_id2efm_ids, id2pattern, r_ids, reversible_r_ids, output_file, sorter=SIMPLE_PATTERN_SORTER):
    """
    Serializes patterns to a file, one pattern per line. Patterns are represented as ids of the active reactions
    (for the reversed reactions, the id is preceded by minus), e.g. -R1 R3 -R7 R11 R25.
    Patterns are sorted according to the sorter function.

    :param p_id2efm_ids: dict, {pattern id: ids of EFMs containing this pattern}

    :param id2pattern: dict, {pattern id: pattern}. Patterns are represented in binary form.

    A binary representation of a pattern is a list of integers whose binary representations
    correspond to the reactions that are active in this pattern: if the reaction is active,
    the corresponding bit is set to 1.
    If the total number of reactions in the model is larger that the number of bits in an int, several ints are used.

    Example: For a model containing reactions r1, r2(reversible), r3(reversible), r4, r5,
    a pattern: r1, r2, r3, r5 would be represented as [77], as the binary representation of 77 is '1001101'
    that corresponds to '1 r5, 0 r4, 0 -r3, 1 r3, 1 -r2, 0 r2, 1 r1', where -ri is the reversed version of ri.

    :param r_ids: ordered collection of reaction ids (strings).

    :param reversible_r_ids: set of ids of reversible reactions (strings).

    :param output_file: path to the file where the patterns should be saved

    :param sorter: a function that given a pattern id return a value that will be used to sort a collection of patterns
    """
    int_size = math.log(sys.maxint) / math.log(2)

    def log_pattern(p_id, f):
        pattern = id2pattern[p_id]
        p_length = get_binary_efm_length(pattern)
        efm_len = len(p_id2efm_ids[p_id])
        p_string = '\t'.join('%s%s' % ('-' if coeff < 0 else '', r_id)
                             for (r_id, coeff) in
                             binary2efm(pattern, r_ids, reversible_r_ids, int_size, binary=True).iteritems())
        f.write("(%d, %d)\t%s\n" % (p_length, efm_len, p_string))

    with open(output_file, 'w+') as f:
        f.write("!(p_length, efm_num)\tpattern\n")
        for p_id in sorted(p_id2efm_ids.iterkeys(), key=sorter):
            log_pattern(p_id, f)


def serialize_pattern_hierarchy(p_id2efm_ids, id2pattern, child2parents, r_id2i, output_file):
    """
    Serializes a pattern hierarchy to a file, one pattern per line.
    Patterns are represented as ids of the active reactions
    (for the reversed reactions, the id is preceded by minus), e.g. -R1 R3 -R7 R11 R25.
    Subpatterns are placed on the lines following their superpattern's line, preceded by a tabulation.

    :param p_id2efm_ids: dict, {pattern id: ids of EFMs containing this pattern}

    :param id2pattern: dict, {pattern id: pattern}. Patterns are represented in binary form.

    A binary representation of a pattern is a list of integers whose binary representations
    correspond to the reactions that are active in this pattern: if the reaction is active,
    the corresponding bit is set to 1.
    If the total number of reactions in the model is larger that the number of bits in an int, several ints are used.

    Example: For a model containing reactions r1, r2(reversible), r3(reversible), r4, r5,
    a pattern: r1, r2, r3, r5 would be represented as [77], as the binary representation of 77 is '1001101'
    that corresponds to '1 r5, 0 r4, 0 -r3, 1 r3, 1 -r2, 0 r2, 1 r1', where -ri is the reversed version of ri.

    :param r_id2i: dict, {reaction_id: positive integer, associated with this reaction}

    :param output_file: path to the file where the patterns should be saved

    :param child2parents: dict, {subpattern_id: set of ids of its superpatterns}
    """
    parent2children = defaultdict(set)
    for ch_p_id, par_p_ids in child2parents.iteritems():
        for par_p_id in par_p_ids:
            parent2children[par_p_id].add(ch_p_id)
    roots = [par_p_id for par_p_id in parent2children.iterkeys() if par_p_id not in child2parents]
    i2r_id = {i: r_id for (r_id, i) in r_id2i.iteritems()}
    i2r_id.update({-i: '-%s' % r_id for (r_id, i) in r_id2i.iteritems()})

    def log_pattern(p_id, level):
        with open(output_file, 'a+') as f:
            f.write("%s%s (%d)\n" % ('\t' * level, " ".join(i2r_id[i] for i in id2pattern[p_id]),
                                     len(p_id2efm_ids[p_id])))

    dfs(roots, parent2children, level=0, process=lambda it, level: log_pattern(it, level))


def dfs(roots, element2children, process, level):
    """
    Preorder depth-first search.
    :param roots: an iterable of root elements
    :param element2children: dict, {element: iterable of child elements}
    :param process: the operation to be performed on each element, lambda element, level: operation(element, level)
    :param level: the level of the current root in the complete tree.
    """
    for it in roots:
        process(it, level)
        if it in element2children:
            dfs(element2children[it], element2children, process, level + 1)