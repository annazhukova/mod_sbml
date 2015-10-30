from mod_sbml.annotation.chebi.chebi_serializer import get_chebi
from mod_sbml.annotation.gene_ontology.go_serializer import get_go
from mod_sbml.annotation.kegg.kegg_annotator import annotate_compounds, annotate_reactions, annotate_pathways
from mod_sbml.annotation.chebi.chebi_annotator import annotate_metabolites
from mod_sbml.annotation.gene_ontology.go_annotator import annotate_compartments
from mod_sbml.onto import parse_simple

__author__ = 'anna'


def annotate(model, compartments=True, metabolites=True, reactions=True, pathways=True, pw_threshold=0.5, org=None,
             chebi=None):
    if compartments:
        go = parse_simple(get_go())
        annotate_compartments(model, go)
    if metabolites or reactions or pathways:
        if not chebi:
            chebi = parse_simple(get_chebi())
        annotate_metabolites(model, chebi)
        annotate_compounds(model, chebi)
    if reactions or pathways:
        annotate_reactions(model)
    if pathways:
        annotate_pathways(model, threshold=pw_threshold, org=org)


