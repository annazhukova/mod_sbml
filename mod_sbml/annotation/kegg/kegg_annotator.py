from collections import defaultdict
import libsbml
from mod_sbml.annotation.kegg.pathway_manager import get_relevant_pathway_info
from mod_sbml.annotation.kegg.reaction_manager import get_compounds2rn
from mod_sbml.sbml.sbml_manager import get_reactants, get_products
from mod_sbml.annotation.rdf_annotation_helper import get_is_annotations, get_annotations, add_annotation

__author__ = 'anna'

KEGG_PROTON = 'C00080'


def get_kegg_r_id(r):
    for annotation in get_is_annotations(r):
        if annotation.find("kegg.reaction") != -1:
            return annotation.replace("kegg.reaction:", '')
    return None


def get_kegg_m_id(m):
    for annotation in get_is_annotations(m):
        if annotation.find("kegg.compound") != -1:
            return annotation.replace("kegg.compound:", '')
    return None


def find_kegg_m_id(m):
    st = m.getSpeciesType()
    kegg = None
    if st:
        kegg = get_kegg_m_id(st)
    if not kegg:
        kegg = get_kegg_m_id(m)
    return kegg


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


def annotate_with_pathways(org, model, kegg_r_ids, kegg_r_id2r_ids, threshold=0.5):
    pw2name_rs_ratio = get_relevant_pathway_info(org, {"rn:" + r_id for r_id in kegg_r_ids}, threshold=threshold)
    for pw, (name, rns, ratio) in pw2name_rs_ratio.iteritems():
        pw = pw.replace("pathway:", "")
        for kegg_r_id in rns:
            kegg_r_id = kegg_r_id.replace("rn:", "")
            for r_id in kegg_r_id2r_ids(kegg_r_id):
                r = model.getElementBySId(r_id)
                if r:
                    add_annotation(r, libsbml.BQB_IS_PART_OF, "http://identifiers.org/kegg.pathway/%s" % pw)
                    r.appendNotes("<html:body><html:p>SUBSYSTEM: %s</html:p></html:body>" % name)


def infer_reaction_kegg_from_compounds_kegg(model):
    rs_ps2rn = get_compounds2rn()
    m_id2kegg = {m.id: get_kegg_m_id(m) for m in model.getListOfSpecies()}
    for r in model.getListOfReactions():
        if get_kegg_r_id(r):
            continue
        rs, ps = [m_id2kegg[m_id] for m_id in get_reactants(r)], [m_id2kegg[m_id] for m_id in get_products(r)]
        rs, ps = tuple(sorted(rs)), tuple(sorted(ps))
        if rs > ps:
            rs, ps = ps, rs
        if (rs, ps) in rs_ps2rn:
            kegg = rs_ps2rn[(rs, ps)]
            add_annotation(r, libsbml.BQB_IS, "http://identifiers.org/kegg.reaction/%s" % kegg)