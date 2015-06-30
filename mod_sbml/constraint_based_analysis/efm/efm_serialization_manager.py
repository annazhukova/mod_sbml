import logging
import os

import libsbml
import openpyxl
from mod_sbml.utils.misc import invert_map

from mod_sbml.constraint_based_analysis.efm.EFM import EFM
from mod_sbml.sbml.sbml_manager import get_r_comps, get_kegg_r_id, reverse_reaction
from mod_sbml.sbml.submodel_manager import submodel
from mod_sbml.serialization.serialization_manager import get_sbml_r_formula
from mod_sbml.serialization.xlsx_helper import BASIC_STYLE, save_data

SIMPLE_PATTERN_SORTER = lambda p_id: -p_id

__author__ = 'anna'

basic_r_style = lambda r_id: BASIC_STYLE

SORT_BY_PATTERN_LENGTH = 1
SORT_BY_EFM_NUMBER = 0
SORT_BY_WEIGHTED_PRODUCT_OF_EFM_NUMBER_AND_PATTERN_LENGTH = 2


def get_pattern_sorter(id2pattern, p_id2efm_ids, sort=SORT_BY_EFM_NUMBER):
    if not id2pattern:
        return SIMPLE_PATTERN_SORTER
    if SORT_BY_PATTERN_LENGTH == sort:
        return lambda p_id: -len(id2pattern[p_id])
    elif SORT_BY_EFM_NUMBER == sort:
        return lambda p_id: -len(p_id2efm_ids[p_id])
    elif SORT_BY_WEIGHTED_PRODUCT_OF_EFM_NUMBER_AND_PATTERN_LENGTH == sort:
        max_pattern_len = max(len(pattern) for pattern in id2pattern.itervalues())
        min_pattern_len = min(len(pattern) for pattern in id2pattern.itervalues())
        max_efm_number = max(len(efm_ids) for efm_ids in p_id2efm_ids.itervalues())
        min_efm_num = min(len(efm_ids) for efm_ids in p_id2efm_ids.itervalues())
        return lambda p_id: \
            -(1.0 * (len(id2pattern[p_id]) - min_pattern_len) / (max_pattern_len - min_pattern_len)
              if max_pattern_len != min_pattern_len else 1) * \
            ((1.0 * len(p_id2efm_ids[p_id]) - min_efm_num) / (max_efm_number - min_efm_num)
             if max_efm_number != min_efm_num else 1)
    else:
        logging.error('Unknown pattern sorting function %s. Will sort by pattern id.' % sort)
        return SIMPLE_PATTERN_SORTER


def serialize_efms_xslx(sbml, id2efm, path, r_id2style=basic_r_style):
    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()
    wb = openpyxl.Workbook()
    i = 0
    for efm_id, efm in sorted(id2efm.iteritems(), key=lambda (_, efm): len(efm)):
        r_id2coefficient = efm.to_r_id2coeff()
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
        save_data(["Coefficient", "Id", "Name", "Formula", "Kegg", "Compartments"], data=data,
                  ws_name="EM%d_%d" % (efm_id, len(efm)), wb=wb, ws_index=i, styles=styles)
        i += 1
    wb.save(path)


def serialize_efms_txt(id2efm, path):
    with open(path, 'w+') as f:
        f.write('Found %d EFMs:\n' % len(id2efm))
        for efm_id in sorted(id2efm.iterkeys()):
            efm = id2efm[efm_id]
            f.write('%d\t%s\n' % (efm_id, efm))


def serialize_n_most_effective_efms_txt(id2efm, id2efficiency, n, path):
    with open(path, 'w+') as f:
        f.write('%d most effective EFMs:\n' % n)
        for efm_id, efficiency in sorted(id2efficiency.iteritems(), key=lambda (_, eff): -eff)[0:n]:
            efm = id2efm[efm_id]
            f.write('Efficiency=%g\tid=%d\t%s\n' % (efficiency, efm_id, efm))


def read_efms(output_efm_file, r_ids, rev_r_ids):
    id2efm = {}
    with open(output_efm_file, 'r') as f:
        for line in f:
            line = line.replace('\n', '').strip()
            # header line
            if line.find('Found ') == 0:
                continue
            efm = line.split('\t')
            if not efm:
                continue
            efm_id = int(efm[0])
            efm = efm[1:]
            id2efm[efm_id] = EFM(r_id2coeff={r_id: float(coeff) for (coeff, r_id) in (it.split(' ') for it in efm)},
                                 r_ids=r_ids, rev_r_ids=rev_r_ids)
    return id2efm


def efm2sbml(id2efm, directory, sbml, process_sbml=lambda sbml_file, path: True,
             limit=None, get_file_path=lambda f: f, all_efms_file=None):
    count = 0
    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()
    for efm_id, efm in sorted(id2efm.iteritems(), key=lambda (efm_id, efm): len(efm)):
        r_id2coeff = efm.to_r_id2coeff()

        def r_updater(r):
            coeff = r_id2coeff[r.id]
            if coeff < 0:
                reverse_reaction(r)
                r.setName(('%g -%s' % (-coeff, r.id)) if coeff != -1 else ('-%s' % r.id))
            else:
                r.setName(('%g %s' % (coeff, r.id)) if coeff != 1 else ('%s' % r.id))

        efm_name = 'EFM_%d_of_len_%d' % (efm_id, len(efm))
        new_sbml = os.path.join(directory, '%s.xml' % efm_name)
        r_ids2sbml(r_id2coeff.keys(), sbml, new_sbml, efm_name, r_updater)

        efm_txt = os.path.join(directory, '%s.txt' % efm_name)
        serialize_efm(efm_id, r_id2coeff, model, efm_txt)

        efm_file_string = \
            '<p>Detailed description of all <b>%d</b> found EFMs is <a href="%s" target="_blank">here</a>.</p>' \
            % (len(id2efm), get_file_path(all_efms_file)) if all_efms_file else ''
        description = \
            '''<h3>Elementary Flux Mode Analysis (EFMA)</h3>
               <p>EFM %d of length <b>%d</b>
               (<a href="%s" target="_blank">detailed description</a>, <a href="%s" download>SBML submodel</a>) is visualised below.</p>
               %s
            ''' % (efm_id, len(efm), get_file_path(efm_txt), get_file_path(new_sbml), efm_file_string)
        process_sbml(new_sbml, efm_name, description)
        count += 1
        if limit and count >= limit:
            break


def pattern2sbml(id2pattern, id2efm_ids, directory, sbml, process_sbml=lambda sbml_file, path: True,
                 limit=None, get_file_path=lambda f: f, all_patterns_file=None, total_efm_num=0):
    count = 0
    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()
    for p_id, pattern in sorted(id2pattern.iteritems(), key=lambda (p_id, _): -len(id2efm_ids[p_id])):
        r_id2coeff = pattern.to_r_id2coeff()

        def r_updater(r):
            coeff = r_id2coeff[r.id]
            if coeff < 0:
                reverse_reaction(r)
                r.setName(('%g -%s' % (-coeff, r.id)) if coeff != -1 else ('-%s' % r.id))
            else:
                r.setName(('%g %s' % (coeff, r.id)) if coeff != 1 else ('%s' % r.id))

        pattern_name = 'Pattern_%d_of_len_%d_of_%d_EFMs' % (p_id, len(pattern), len(id2efm_ids[p_id]))
        new_sbml = os.path.join(directory, '%s.xml' % pattern_name)
        r_ids2sbml(r_id2coeff.keys(), sbml, new_sbml, pattern_name, r_updater)
        pattern_txt = os.path.join(directory, '%s.txt' % pattern_name)
        serialize_pattern(p_id, r_id2coeff, id2efm_ids[p_id], model, pattern_txt)
        p_file_string = '<p>Detailed description of all found patterns is <a href="%s" target="_blank">here</a>.</p>' \
                        % get_file_path(all_patterns_file) if all_patterns_file else ''
        description = \
            '''<h3>Elementary Flux Mode Analysis (EFMA): common patterns</h3>
               <p>Pattern %d of length <b>%d</b> common to <b>%d</b> (out of <b>%d</b>) EFMs
               (<a href="%s" target="_blank">detailed description</a>, <a href="%s" download>SBML submodel</a>) is visualised below.</p>
               %s
            ''' % (p_id, len(pattern), len(id2efm_ids[p_id]), total_efm_num, get_file_path(pattern_txt),
                   get_file_path(new_sbml), p_file_string)
        process_sbml(new_sbml, pattern_name, description)
        count += 1
        if limit and count >= limit:
            break


def r_ids2sbml(r_ids, sbml, out_sbml, suffix='', r_updater=lambda r: True):
    doc = libsbml.SBMLReader().readSBML(sbml)
    model = doc.getModel()
    submodel(r_ids, model)
    model.setId('%s_%s' % (model.getId(), suffix))
    model.setName('%s_%s' % (model.getName(), suffix))
    for r_id in r_ids:
        r = model.getReaction(r_id)
        if r:
            r_updater(r)
    libsbml.SBMLWriter().writeSBMLToFile(doc, out_sbml)


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


def serialize_important_reactions(r_id2efm_ids, model, output_file, imp_rn_threshold=0):
    important_r_ids = {r_id[1:] if '-' == r_id[0] else r_id for r_id in r_id2efm_ids.iterkeys()}
    with open(output_file, 'w+') as f:
        f.write('Found %d reactions (including %d unique) that participate in more than %d EFMs\n'
                % (len(r_id2efm_ids), len(important_r_ids), imp_rn_threshold))
        for r_id, efm_ids in sorted(r_id2efm_ids.iteritems(), key=lambda (_, efm_ids): -len(efm_ids)):
            f.write('(%d)\t%s\t%s\n' % (len(efm_ids), r_id,
                                    get_sbml_r_formula(model, model.getReaction(r_id[1:] if '-' == r_id[0] else r_id))))
    return important_r_ids


def serialize_pattern(p_id, r_id2coeff, efm_ids, model, output_file):
    """
    Serializes a pattern to a file.

    :param p_id: pattern id

    :param r_id2coeff: pattern, represented as a dictionary {reaction id: coefficient}.

    :param output_file: path to the file where the pattern should be saved

    :param efm_ids: a collection of ids of EFMs that contain this pattern
    """
    with open(output_file, 'w+') as f:
        f.write('Pattern %d of length %d common for %d following EFMs:\n%s\n------------------\n'
                % (p_id, len(r_id2coeff), len(efm_ids), ', '.join((str(it) for it in efm_ids))))
        coeff2r_id = invert_map(r_id2coeff)
        for coeff, r_ids in sorted(coeff2r_id.iteritems(), key=lambda (coeff, _): (-abs(coeff), -coeff)):
            for r_id in sorted(r_ids):
                f.write('%g\t%s:\t%s\n' % (coeff, r_id, get_sbml_r_formula(model, model.getReaction(r_id), comp=False, id=True)))
            f.write('\n')


def serialize_efm(efm_id, r_id2coeff, model, output_file):
    """
    Serializes an EFM to a file.

    :param efm_id: EFM id

    :param r_id2coeff: EFM, represented as a dictionary {reaction id: coefficient}.

    :param output_file: path to the file where the EFM should be saved
    """
    with open(output_file, 'w+') as f:
        f.write('EFM %d of length %d\n------------------\n' % (efm_id, len(r_id2coeff)))
        coeff2r_id = invert_map(r_id2coeff)
        for coeff, r_ids in sorted(coeff2r_id.iteritems(), key=lambda (coeff, _): (-abs(coeff), -coeff)):
            for r_id in sorted(r_ids):
                f.write('%g\t%s:\t%s\n' % (coeff, r_id, get_sbml_r_formula(model, model.getReaction(r_id), comp=False, id=True)))
            f.write('\n')

def serialize_patterns(p_id2efm_ids, id2pattern, output_file, min_pattern_len, sorter=SIMPLE_PATTERN_SORTER):
    """
    Serializes patterns to a file, one pattern per line. Patterns are represented as ids of the active reactions
    (for the reversed reactions, the id is preceded by minus), e.g. -R1 R3 -R7 R11 R25.
    Patterns are sorted according to the sorter function.

    :param p_id2efm_ids: dict, {pattern id: ids of EFMs containing this pattern}

    :param id2pattern: dict, {pattern id: pattern}.

    :param output_file: path to the file where the patterns should be saved

    :param sorter: a function that given a pattern id return a value that will be used to sort a collection of patterns
    """

    def log_pattern(p_id, f):
        pattern = id2pattern[p_id]
        p_length = len(pattern)
        efm_len = len(p_id2efm_ids[p_id])
        p_string = pattern.to_string(binary=True)
        f.write("id:%d\t(len=%d, num of EFMs=%d)\t%s\n" % (p_id, p_length, efm_len, p_string))
        f.write("\tEFMs:\t%s\n" % ', '.join(str(efm_id) for efm_id in p_id2efm_ids[p_id]))

    with open(output_file, 'w+') as f:
        f.write('Found %d patterns of len >= %d.\n'
                % (len(id2pattern), min_pattern_len))
        for p_id in sorted(p_id2efm_ids.iterkeys(), key=sorter):
            log_pattern(p_id, f)
