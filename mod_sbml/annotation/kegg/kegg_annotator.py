from collections import defaultdict
from itertools import chain
import logging

import libsbml

from mod_sbml.annotation.chebi.chebi_annotator import get_chebi_id
from mod_sbml.annotation.kegg.pathway_manager import get_relevant_pathway_info
from mod_sbml.annotation.kegg.reaction_manager import get_compounds2rn, get_kegg_r_id_by_kegg_m_ids
from mod_sbml.sbml.sbml_manager import get_reactants, get_products, get_subsystem2r_ids
from mod_sbml.annotation.rdf_annotation_helper import get_is_annotations, add_annotation, get_annotations

KEGG_REACTION_PREFIX = "kegg.reaction"

KEGG_COMPOUND_PREFIX = "kegg.compound"
KEGG_PATHWAY_PREFIX = "kegg.pathway"

__author__ = 'anna'

KEGG_PROTON = 'C00080'


def get_kegg_r_id(r):
    for annotation in get_is_annotations(r):
        if annotation.find(KEGG_REACTION_PREFIX) != -1:
            return annotation.replace("%s:" % KEGG_REACTION_PREFIX, '')
    return None


def get_kegg_m_id(m):
    for annotation in get_is_annotations(m):
        if annotation.find(KEGG_COMPOUND_PREFIX) != -1:
            return annotation.replace("%s:" % KEGG_COMPOUND_PREFIX, '')
    return None


def infer_kegg_m_id(m, chebi):
    if isinstance(m, libsbml.SpeciesType):
        st = m.getSpeciesType()
        if st:
            kegg_id = get_kegg_m_id(st)
            if kegg_id:
                return kegg_id
    chebi_id = get_chebi_id(m)
    if not chebi_id or not chebi:
        return None
    term = chebi.get_term(chebi_id)
    if term:
        kegg_ids = sorted(term.get_kegg_ids())
        if kegg_ids:
            return kegg_ids[0].upper()


def annotate_compounds(model, chebi=None):
    for m in model.getListOfSpecies():
        if get_kegg_m_id(m):
            continue
        kegg_id = infer_kegg_m_id(m, chebi)
        if kegg_id:
            add_annotation(m, libsbml.BQB_IS, kegg_id, KEGG_COMPOUND_PREFIX)


def get_pathway2r_ids(sbml=None, model=None):
    pw2r_ids = defaultdict(set)
    no_pw_r_ids = set()
    if not model and not sbml:
        raise ValueError("Either sbml or model parameter should be specified")
    if not model:
        input_doc = libsbml.SBMLReader().readSBML(sbml)
        model = input_doc.getModel()
    for r in model.getListOfReactions():
        found = False
        for annotation in get_annotations(r, libsbml.BQB_IS_PART_OF):
            if annotation.find("path:") != -1:
                pw2r_ids[annotation.replace("path:", '')].add(r.getId())
                found = True
        if not found:
            no_pw_r_ids.add(r.getId())
    return pw2r_ids, no_pw_r_ids


def get_kegg_r_id2r_ids(model):
    kegg_r_id2r_ids = defaultdict(set)
    for r in model.getListOfReactions():
        k_r_id = get_kegg_r_id(r)
        if k_r_id:
            kegg_r_id2r_ids[k_r_id].add(r.getId())
    return kegg_r_id2r_ids


def get_kegg_m_id2m_ids(model):
    kegg_m_id2m_ids = defaultdict(set)
    for m in model.getListOfSpecies():
        k_m_id = get_kegg_m_id(m)
        if k_m_id:
            kegg_m_id2m_ids[k_m_id].add(m.getId())
    return kegg_m_id2m_ids


def annotate_pathways(model, threshold=0.5, org='map'):
    if get_subsystem2r_ids(model=model)[0]:
        return
    if not org:
        org = 'map'
    kegg_r_id2r_ids = get_kegg_r_id2r_ids(model)
    kegg_r_ids = set(kegg_r_id2r_ids.keys())
    pw2name_rs_ratio = []
    try:
        pw2name_rs_ratio = get_relevant_pathway_info(org, {"rn:" + r_id for r_id in kegg_r_ids}, threshold=threshold)
    except Exception as e:
        logging.error('Did not manage to infer pathways due to %s' % e.message)
    for pw, (name, rns, ratio) in pw2name_rs_ratio.items():
        pw = pw.replace("pathway:", "")
        for kegg_r_id in rns:
            kegg_r_id = kegg_r_id.replace("rn:", "")
            for r_id in kegg_r_id2r_ids[kegg_r_id]:
                r = model.getElementBySId(r_id)
                if r:
                    add_annotation(r, libsbml.BQB_IS_PART_OF, pw, KEGG_PATHWAY_PREFIX)
                    r.appendNotes("<html:body><html:p>SUBSYSTEM: %s</html:p></html:body>" % name)


def annotate_reactions(model):
    m_id2kegg = None
    for r in model.getListOfReactions():
        if get_kegg_r_id(r):
            continue
        if m_id2kegg is None:
            m_id2kegg = {m.getId(): get_kegg_m_id(m) for m in model.getListOfSpecies()}
        kegg_id = infer_r_kegg_id(r, m_id2kegg)
        if kegg_id:
            add_annotation(r, libsbml.BQB_IS, kegg_id, KEGG_REACTION_PREFIX)


def infer_r_kegg_id(r, m_id2kegg):
    ms = {m_id2kegg[m_id] for m_id in chain(get_reactants(r), get_products(r))}
    if None in ms:
        return None
    try:
        kegg_ids = get_kegg_r_id_by_kegg_m_ids(ms)
        if kegg_ids and len(kegg_ids) == 1:
            return kegg_ids.pop()
    except Exception as e:
        logging.error(e.message)
    return None
