import libsbml

from mod_sbml.onto.chebi_annotator import get_chebi, get_species_to_chebi
from mod_sbml.kegg.pathway_manager import get_relevant_pathway_info
from mod_sbml.kegg.reaction_manager import get_compounds2rn
from model_comparison.model_matcher import match_ms, match_rs
from mod_sbml.onto.obo_ontology import parse
from mod_sbml.sbml.sbml_manager import get_kegg_r_id, add_annotation, get_kegg_m_id, get_reactants, get_products


__author__ = 'anna'


def copy_m_kegg_annotations(model, template_model):
    chebi = parse(get_chebi())
    m_id2chebi = get_species_to_chebi(model, chebi)
    t_m_id2chebi = get_species_to_chebi(template_model, chebi)
    m_id2group, t_m_id2group, m_i2i, m_group2info, t_m_group2info, group2m_ids, _ = \
        match_ms(model, template_model, m_id2chebi, t_m_id2chebi, chebi)

    for i, m_ids in group2m_ids.iteritems():
        kegg_ids, _, _, _ = m_group2info[i]
        if not kegg_ids and i in m_i2i:
            kegg_ids, _, _, _ = t_m_group2info[m_i2i[i]]
        if not kegg_ids:
            continue
        kegg = list(kegg_ids)[0].upper()
        for m_id in m_ids:
            m = model.getSpecies(m_id)
            old_kegg = get_kegg_m_id(m)
            if not old_kegg:
                add_annotation(m, libsbml.BQB_IS, "http://identifiers.org/kegg.compound/%s" % kegg)

    return m_i2i, m_id2group, t_m_id2group


def copy_kegg_annotations(model, template_model):
    m_i2i, m_id2group, t_m_id2group = copy_m_kegg_annotations(model, template_model)

    r_id2group, t_r_id2group, r_i2i, r_group2info, t_r_group2info, group2r_ids, _ = \
        match_rs(model, template_model, m_id2group, t_m_id2group, m_i2i)

    for i, r_ids in group2r_ids.iteritems():
        kegg_ids, _ = r_group2info[i]
        if not kegg_ids and i in r_i2i:
            kegg_ids, _ = t_r_group2info[r_i2i[i]]
        if not kegg_ids:
            continue
        kegg = list(kegg_ids)[0].upper()
        for r_id in r_ids:
            r = model.getReaction(r_id)
            old_kegg = get_kegg_r_id(r)
            if not old_kegg:
                add_annotation(r, libsbml.BQB_IS, "http://identifiers.org/kegg.reaction/%s" % kegg)


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