from collections import defaultdict

# from mod_sbml.annotation.chebi.chebi_annotator import get_species_to_chebi
from mod_sbml.annotation.kegg.pathway_manager import get_pw_name
from mod_sbml.annotation.kegg.kegg_annotator import get_pathway2r_ids
# from mod_sbml.onto.onto_getter import get_chebi
# from mod_sbml.annotation.kegg.reaction_manager import get_formula2kegg_compound
# from model_comparison.model_matcher import group_metabolites, group_rs
# from mod_sbml.onto.obo_ontology import parse
from mod_sbml.sbml.sbml_manager import get_metabolites, get_subsystem2r_ids, \
    get_reactants, get_products


__author__ = 'anna'


def model_statistics(model, org, groups=True, pathways=True, extracellular=None):
    print("Compartments: ", model.getNumCompartments(), "Metabolites: ", model.getNumSpecies(), "Reactions: ",
          model.getNumReactions())
    print("Compartments: %s" % ", ".join([c.name for c in model.getListOfCompartments()]))
    c2reactions = {}
    s_id2rs = defaultdict(list)
    s_id2r_rs = defaultdict(list)
    s_id2p_rs = defaultdict(list)
    transport = []
    for r in model.getListOfReactions():
        rs, ps = get_reactants(r), get_products(r)
        ms = get_metabolites(r)
        for m in rs:
            s_id2r_rs[m].append(r.id)
        for m in ps:
            s_id2p_rs[m].append(r.id)
        for m in ms:
            s_id2rs[m].append(r)
        c_ids = tuple(sorted({model.getSpecies(s_id).getCompartment() for s_id in ms}))
        if len(c_ids) > 1:
            transport.append(r.id)
        if c_ids not in c2reactions:
            c2reactions[c_ids] = 1
        else:
            c2reactions[c_ids] += 1

    boundary_ms = {m for m in model.getListOfSpecies() if m.id in s_id2rs and m.getBoundaryCondition()}
    if not boundary_ms and extracellular:
        for m in (m for m in model.getListOfSpecies() if m.getCompartment() == extracellular):
            for r in s_id2rs[m.id]:
                if not set(get_reactants(r)) or not set(get_products(r)):
                    boundary_ms.add(m)
                    break
    print("Boundary metabolites: ", len(boundary_ms))
    boundary_ms_ids = {s.id for s in boundary_ms}
    blocked_ms = {s_id for s_id in s_id2r_rs if s_id not in s_id2p_rs and s_id not in boundary_ms_ids} | \
                 {s_id for s_id in s_id2p_rs if s_id not in s_id2r_rs and s_id not in boundary_ms_ids}
    print("Real metabolites: ", len(s_id2rs.keys()))
    print("Blocked metabolites: ", len(blocked_ms))
    print("Transport reactions: ", len(transport))
    print("-------------Reaction distribution-----------------")
    for c_ids in sorted(c2reactions.iterkeys()):
        print([model.getCompartment(c_id).name for c_id in c_ids], c2reactions[c_ids])
    print("---------------------------------------------------")

    # # Element groups
    # if groups:
    #     chebi = parse(get_chebi())
    #     m_id2chebi = get_species_to_chebi(model, chebi)
    #     formula2kegg = get_formula2kegg_compound()
    #     group2m_ids, m_id2group, m_group2info = group_metabolites(model, m_id2chebi, chebi, formula2kegg)
    #     print("Metabolite types: ", len(m_group2info.keys()))
    #     group2r_ids, r_id2group, r_group2info = group_rs(model, lambda m_id: m_id2group[m_id])
    #     print("Reaction types: ", len(r_group2info.keys()))
    #     transport_groups = {r_id2group[r_id] for r_id in transport}
    #     print("Transport types: ", len(transport_groups))

    # Pathways
    if pathways:
        pw2r_ids, o_r_ids = get_pathway2r_ids(model=model)
        if not pw2r_ids.keys():
            pw2r_ids, o_r_ids = get_subsystem2r_ids(model=model)
        pw_names = []
        for pw in pw2r_ids.keys():
            pw_name = get_pw_name(org, pw)
            if pw_name:
                pw = "%s (%s)" % (pw_name, pw)
            pw_names.append(pw)
        print ("Pathways", sorted(pw_names))

