import libsbml
from annotation.annotation_manager import infer_reaction_kegg_from_compounds_kegg
from chebi_annotator import get_chebi, get_species_to_chebi
from kegg.reaction_manager import get_formula2kegg_compound
from model_comparison.model_matcher import group_metabolites
from obo_ontology import parse
from sbml_manager import get_kegg_m_id, add_annotation

__author__ = 'anna'


def annotate_model_with_kegg(model):
    chebi = parse(get_chebi())
    m_id2chebi = get_species_to_chebi(model, chebi)
    formula2kegg = get_formula2kegg_compound()
    group2m_ids, m_id2group, group2info = group_metabolites(model, m_id2chebi, chebi, formula2kegg)
    for gr, m_ids in group2m_ids.iteritems():
        kegg_ids, _, _, _ = group2info[gr]
        if not kegg_ids:
            continue
        kegg = list(kegg_ids)[0].upper()
        for m_id in m_ids:
            m = model.getSpecies(m_id)
            old_kegg = get_kegg_m_id(m)
            if not old_kegg:
                add_annotation(m, libsbml.BQB_IS, "http://identifiers.org/kegg.compound/%s" % kegg)

    infer_reaction_kegg_from_compounds_kegg(model)
