from cobra.io.sbml import write_cobra_model_to_sbml_file

__author__ = 'anna'


def fix_names(model):
    for species in model.metabolites:
        species.name = species.name.replace(species.id, '').strip()


def keep_only_reactions(model, r_ids, file_path):
    to_delete = [r for r in model.reactions if r.id not in r_ids]
    for r in to_delete:
        model.reactions.remove(r.id)
    write_cobra_model_to_sbml_file(model, file_path, 2, 4)


def get_reaction_ids_by_participant_metabolite_id(model, s_id):
    return {r.id for r in model.reactions if
            s_id in {m.id for m in r.metabolites}}


def get_reaction_ids_between_metabolite_ids(model, s_id, s_other_id):
    result = set()
    for r in model.reactions:
        rs, ps = {m.id for m in r.products}, {m.id for m in r.reactants}
        if s_id in rs and s_other_id in ps or s_id in ps and s_other_id in rs:
            result.add(r.id)
    return result


def get_reaction_ids_between_metabolite_id_sets(model, s_ids, s_other_ids):
    result = set()
    for r in model.reactions:
        rs, ps = {m.id for m in r.products}, {m.id for m in r.reactants}
        if s_ids & rs and s_other_ids & ps or s_ids & ps and s_other_ids & rs:
            result.add(r.id)
    return result


def get_transport_reactions(model, c_id=None):
    result = set()
    for r in model.reactions:
        c_ids = set(r.get_compartments())
        if len(c_ids) > 1 and (c_id is None or c_id in c_ids):
            result.add(r.id)
    return result


def get_boundary_reactions(model):
    result = set()
    for r in model.reactions:
        if len(r.reactants) == 0 or len(r.products) == 0:
            result.add(r.id)
    return result


def format_r_id(r_id, remove_prefix=True):
    if remove_prefix and r_id.startswith('R_'):
        r_id = r_id[2:]
    if not remove_prefix and not r_id.startswith('R_'):
        return 'R_' + r_id
    return r_id
