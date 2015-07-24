from collections import defaultdict
import libsbml
# from kegg.reaction_manager import get_formula2kegg_compound
#
# from mod_sbml.annotation.chebi.chebi_annotator import get_species_to_chebi
# from mod_sbml.onto.onto_getter import get_chebi
from mod_sbml.annotation.kegg.pathway_manager import get_relevant_pathway_info
from mod_sbml.annotation.kegg.reaction_manager import get_compounds2rn
# from model_comparison.model_matcher import match_ms, match_rs, group_metabolites
# from mod_sbml.onto.obo_ontology import parse
from mod_sbml.sbml.sbml_manager import get_reactants, get_products
from rdf_annotation_helper import get_is_annotations, get_annotations, add_annotation

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
            if annotation.find("kegg.pathway") != -1:
                pw2r_ids[annotation.replace("kegg.pathway:", '')].add(r.getId())
                found = True
        if not found:
            no_pw_r_ids.add(r.getId())
    # Should find the following pathways:
    # 1. citrate cycle (TCA cycle) (hsa00020);
    # 2. fatty acid metabolism (hsa01212) and 3 sub-pathways:
    #   (a) fatty acid biosynthesis (hsa00061);
    #   (b) fatty acid elongation (hsa00062);
    #   (c) fatty acid degradation (hsa00071);
    # 3. valine, leucine and isoleucine degradation (hsa00280);
    # 4. synthesis and degradation of ketone bodies (hsa00072).

    # Let's remove the sub-pathways
    if "path:hsa01212" in pw2r_ids:
        for key in ("path:hsa00061", "path:hsa00062", "path:hsa00071"):
            if key in pw2r_ids:
                del pw2r_ids[key]

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


# def annotate_model_with_kegg(model):
#     chebi = parse(get_chebi())
#     m_id2chebi = get_species_to_chebi(model, chebi)
#     formula2kegg = get_formula2kegg_compound()
#     group2m_ids, m_id2group, group2info = group_metabolites(model, m_id2chebi, chebi, formula2kegg)
#     for gr, m_ids in group2m_ids.iteritems():
#         kegg_ids, _, _, _ = group2info[gr]
#         if not kegg_ids:
#             continue
#         kegg = list(kegg_ids)[0].upper()
#         for m_id in m_ids:
#             m = model.getSpecies(m_id)
#             old_kegg = get_kegg_m_id(m)
#             if not old_kegg:
#                 add_annotation(m, libsbml.BQB_IS, "http://identifiers.org/kegg.compound/%s" % kegg)
#
#     infer_reaction_kegg_from_compounds_kegg(model)


# def copy_m_kegg_annotations(model, template_model):
#     chebi = parse(get_chebi())
#     m_id2chebi = get_species_to_chebi(model, chebi)
#     t_m_id2chebi = get_species_to_chebi(template_model, chebi)
#     m_id2group, t_m_id2group, m_i2i, m_group2info, t_m_group2info, group2m_ids, _ = \
#         match_ms(model, template_model, m_id2chebi, t_m_id2chebi, chebi)
#
#     for i, m_ids in group2m_ids.iteritems():
#         kegg_ids, _, _, _ = m_group2info[i]
#         if not kegg_ids and i in m_i2i:
#             kegg_ids, _, _, _ = t_m_group2info[m_i2i[i]]
#         if not kegg_ids:
#             continue
#         kegg = list(kegg_ids)[0].upper()
#         for m_id in m_ids:
#             m = model.getSpecies(m_id)
#             old_kegg = get_kegg_m_id(m)
#             if not old_kegg:
#                 add_annotation(m, libsbml.BQB_IS, "http://identifiers.org/kegg.compound/%s" % kegg)
#
#     return m_i2i, m_id2group, t_m_id2group


# def copy_kegg_annotations(model, template_model):
#     m_i2i, m_id2group, t_m_id2group = copy_m_kegg_annotations(model, template_model)
#
#     r_id2group, t_r_id2group, r_i2i, r_group2info, t_r_group2info, group2r_ids, _ = \
#         match_rs(model, template_model, m_id2group, t_m_id2group, m_i2i)
#
#     for i, r_ids in group2r_ids.iteritems():
#         kegg_ids, _ = r_group2info[i]
#         if not kegg_ids and i in r_i2i:
#             kegg_ids, _ = t_r_group2info[r_i2i[i]]
#         if not kegg_ids:
#             continue
#         kegg = list(kegg_ids)[0].upper()
#         for r_id in r_ids:
#             r = model.getReaction(r_id)
#             old_kegg = get_kegg_r_id(r)
#             if not old_kegg:
#                 add_annotation(r, libsbml.BQB_IS, "http://identifiers.org/kegg.reaction/%s" % kegg)


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
                    r.appendNotes("<html:body><html:p>SUBSYSTEM: {1} ({0})</html:p></html:body>".format(pw, name))


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