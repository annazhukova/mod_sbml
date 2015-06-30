from collections import defaultdict
import logging

from cobra.flux_analysis import flux_variability_analysis
from cobra import Reaction
from cobra.flux_analysis.parsimonious import optimize_minimal_flux

from mod_sbml.constraint_based_analysis.cobra_constraint_based_analysis.model_manager import get_transport_reactions, \
    get_boundary_reactions
from mod_sbml.serialization.plot_manager import create_plot, save_fig, create_subplot, initialise_fig
from mod_sbml.serialization.serialization_manager import get_cobra_r_formula

__author__ = 'anna'

ROUND_COEFFICIENT = 6
NUM_STEPS = 20.0


def set_bounds(reaction, value, max_up_b, max_l_b):
    value = abs(value)
    reaction.lower_bound = max(0 - value, 0 - abs(max_l_b))
    reaction.upper_bound = min(value, abs(max_up_b))


def optimise_biomass(model, bm_id, objective_sense='maximize', level=logging.INFO):
    reaction = model.reactions.get_by_id(bm_id)
    model.change_objective([reaction])
    optimize_minimal_flux(model, objective_sense=objective_sense)
    if "optimal" != model.solution.status:
        logging.error("A problem occurred while calculating the solution for %s: %s"
                      % (format_r_name(reaction), model.solution.status))
        return None
    else:
        logging.log(level, "%s: %.4g" % (format_r_name(reaction), model.solution.f))
        return model.solution.f


def round_value(value):
    return round(value, ROUND_COEFFICIENT)


def get_biomass_dependency(model, r_id, r_ids, bm_id, max_bound=None, remove_zeros=False, constrained=False,
                           minimized=False):
    reaction = model.reactions.get_by_id(r_id)
    r_up_b, r_l_b = reaction.upper_bound, reaction.lower_bound
    max_bound = abs(max_bound) if max_bound is not None else max(abs(r_up_b), abs(r_l_b))
    r_id2ys = defaultdict(list)
    step = max_bound / NUM_STEPS
    xs = []

    if constrained:
        r_id2bounds = get_possible_reaction_bounds(model, bm_id, r_ids - {r_id, bm_id}, reaction, xs)
        old_r_id2bounds = constraint_reactions(model, r_id2bounds)

    model.change_objective([model.reactions.get_by_id(bm_id)])

    for x in [i * step for i in range(0, int(NUM_STEPS))]:
        # reaction.lower_bound = -x
        # reaction.upper_bound = -x
        set_bounds(reaction, x, abs(r_l_b), abs(r_up_b))
        model.optimize()
        if model.solution.f:
            xs.append(x)

            if minimized:
                r_id2values = flux_variability_analysis(model, the_reactions=list(set(r_ids) - {r_id, bm_id}))
                r_id_by_flux = [rn_id for (rn_id, val) in r_id2values.iteritems() if
                                rn_id not in {r_id, bm_id} and round_value(val['minimum']) != round_value(
                                    val['maximum'])]
                r_id2bounds = get_reaction_bounds(model, r_id_by_flux)
                while r_id_by_flux:
                    max_flux_r_id = max(r_id_by_flux, key=lambda rn_id: -abs(model.solution.x_dict[rn_id]))
                    min_max = r_id2values[max_flux_r_id]
                    min_v, max_v = min_max['minimum'], min_max['maximum']
                    max_flux_r = model.reactions.get_by_id(max_flux_r_id)
                    flux_value = 0 if min_v * max_v < 0 else min_v if abs(min_v) < abs(max_v) else max_v
                    max_flux_r.lower_bound = flux_value
                    max_flux_r.upper_bound = flux_value
                    model.optimize()
                    r_id2values = flux_variability_analysis(model, the_reactions=r_id_by_flux)
                    r_id_by_flux = [rn_id for (rn_id, val) in r_id2values.iteritems() if
                                    round_value(val['minimum']) != round_value(val['maximum'])]
                constraint_reactions(model, r_id2bounds)

            for r_id in r_ids:
                r_id2ys[r_id].append(round_value(model.solution.x_dict[r_id]))
    reaction.upper_bound = r_up_b
    reaction.lower_bound = r_l_b

    if constrained:
        constraint_reactions(model, old_r_id2bounds)

    if remove_zeros:
        return xs, {r_id: ys for (r_id, ys) in r_id2ys.iteritems() if ys and (len(set(ys)) > 1 or ys[0] != 0)}
    else:
        return xs, r_id2ys


def get_possible_reaction_bounds(model, bm_id, r_ids, r_to_vary, steps):
    model.change_objective([model.reactions.get_by_id(bm_id)])
    model.optimize()
    r_id2min_max = flux_variability_analysis(model, allow_loops=False)
    r_id2bounds = {r_id: (r_id2min_max[r_id]['minimum'], r_id2min_max[r_id]['maximum']) for r_id in r_ids}
    up_b, l_b = abs(r_to_vary.upper_bound), abs(r_to_vary.lower_bound)
    for x in steps:
        set_bounds(r_to_vary, x, up_b, l_b)
        model.optimize()
        variable_r_ids = list(r_id2bounds.keys())
        for r_id in variable_r_ids:
            v_min, v_max = r_id2bounds[r_id]
            r_id2min_max = flux_variability_analysis(model, allow_loops=False, the_reactions=variable_r_ids)[r_id]
            new_v_min, new_v_max = max(r_id2min_max['minimum'], v_min), min(r_id2min_max['maximum'], v_max)
            if new_v_min > new_v_max:
                del r_id2bounds[r_id]
            else:
                r_id2bounds[r_id] = new_v_min, new_v_max
    return r_id2bounds


def constraint_reactions(model, r_id2bounds):
    old_r_id2bounds = {}
    for r_id, (v_min, v_max) in r_id2bounds.iteritems():
        rn = model.reactions.get_by_id(r_id)
        old_r_id2bounds[r_id] = (rn.lower_bound, rn.upper_bound)
        rn.lower_bound = v_min
        rn.upper_bound = v_max
    return old_r_id2bounds


def get_reaction_bounds(model, r_ids):
    r_id2bounds = {}
    for r_id in r_ids:
        r = model.reactions.get_by_id(r_id)
        r_id2bounds[r_id] = r.lower_bound, r.upper_bound
    return r_id2bounds


def get_reaction_fluxes(model, r_ids, bm_id):
    model.change_objective([model.reactions.get_by_id(bm_id)])
    model.optimize()
    r_id2value = {r_id: round_value(model.solution.x_dict[r_id]) for r_id in r_ids}
    r_id2value[bm_id] = round_value(model.solution.f)
    return r_id2value


def get_reactions_with_fluxes_of_value(model, flux_value=0):
    return [r.id for r in model.reactions if round_value(model.solution.x_dict[r.id]) == flux_value]


def get_reactions_with_large_fluxes(model, threshold):
    return [r_id for (r_id, value) in model.solution.x_dict.iteritems() if abs(round_value(value)) > threshold]


def print_fluxes_larger_than_threshold_by_pathway(model, pw2r_ids, o_r_ids, pw2name, threshold=0):
    for pw in sorted(pw2r_ids.iterkeys()):
        r_ids = pw2r_ids[pw]
        logging.info("\n---------------")
        logging.info("%s (%d out of %d)" % (
            pw2name(pw), len([r_id for r_id in r_ids if r_id in model.solution.x_dict and round_value(model.solution.x_dict[r_id])]), len(r_ids)))
        logging.info("---------------\n")
        print_fluxes_larger_than_threshold(model, r_ids=r_ids, threshold=threshold)
    if o_r_ids:
        i_r_ids = get_boundary_reactions(model) & o_r_ids
        t_r_ids = get_transport_reactions(model) & o_r_ids
        ti_r_ids = i_r_ids | t_r_ids
        logging.info("\n---------------")
        logging.info("No pathway assigned: input and transport (%d out of %d)" %
              (len([r_id for r_id in ti_r_ids if r_id in model.solution.x_dict and round_value(model.solution.x_dict[r_id])]), len(ti_r_ids)))
        logging.info("---------------\n")
        print_fluxes_larger_than_threshold(model, r_ids=ti_r_ids, threshold=threshold)
        o_r_ids = o_r_ids - ti_r_ids
        logging.info("\n---------------")
        logging.info("No pathway assigned: other (%d out of %d)" %
              (len([r_id for r_id in o_r_ids if r_id in model.solution.x_dict and round_value(model.solution.x_dict[r_id])]), len(o_r_ids)))
        logging.info("---------------\n")
        print_fluxes_larger_than_threshold(model, r_ids=o_r_ids, threshold=threshold)


def print_fluxes_larger_than_threshold(model, threshold=0, r_ids=None, highlighted_r_ids=None):
    r_id2val = {}
    value2rn = defaultdict(list)
    for r_id, value in model.solution.x_dict.iteritems():
        value = round_value(value)
        if threshold is None or abs(value) > threshold:
            rn = model.reactions.get_by_id(r_id)
            # value2rn[value].append(format_r_name(rn) + ": " + rn.build_reaction_string(True))
            if not r_ids or r_id in r_ids:
                value2rn[value].append(rn.id + ": " + get_cobra_r_formula(rn, comp=False)
                                       + (" *" if highlighted_r_ids and rn.id in highlighted_r_ids else ''))
            r_id2val[r_id] = value

    prev_value = None
    for value in sorted(value2rn.iterkeys(), key=lambda v: abs(v)):
        for rn in sorted(value2rn[value]):
            if prev_value is not None and abs(prev_value) != abs(value):
                logging.info('')
            prev_value = value
            logging.info("%.4g\t%s" % (value, rn))

    return r_id2val


def get_fluxes_larger_than_threshold(model, threshold=0):
    r_id2val = {}
    for r_id, value in model.solution.x_dict.iteritems():
        value = round_value(value)
        if threshold is None or abs(value) > threshold:
            r_id2val[r_id] = value
    return r_id2val


def serialize_fluxes(model, r_id2val, path, r_ids=None, title=None):
    value2rn = defaultdict(list)
    for r_id, value in r_id2val.iteritems():
        if not r_ids or r_id in r_ids:
            rn = model.reactions.get_by_id(r_id)
            value2rn[value].append(r_id + ": " + get_cobra_r_formula(rn, comp=True))

    prev_value = None
    with open(path, 'w+') as f:
        if title:
            f.write(title + '\n')
        for value in sorted(value2rn.iterkeys(), key=lambda v: abs(v)):
            for rn in sorted(value2rn[value]):
                if prev_value is not None and abs(prev_value) != abs(value):
                    f.write('\n')
                prev_value = value
                f.write("%.4g\t%s\n" % (value, rn))


def print_fva(model, rs=None, threshold=None, r_ids=None):
    r_id2min_max = flux_variability_analysis(model, reaction_list=rs)
    values2r_ids = defaultdict(set)
    r_id2bounds = {}
    for r_id, values in r_id2min_max.iteritems():
        min_v, max_v = round_value(values['minimum']), round_value(values['maximum'])
        if threshold is None or abs(min_v) > threshold or abs(max_v) > threshold:
            values2r_ids[(min_v, max_v)].add(r_id)
            r_id2bounds[r_id] = min_v, max_v
    keys = sorted(values2r_ids.iterkeys(), key=lambda (min_v, max_v): (-abs(max_v - min_v), min_v, max_v))
    ess_count, var_count = 0, 0
    logging.info("==============================\nESSENTIAL FLUXES\n==============================")
    for (min_v, max_v) in keys:
        if min_v * max_v > 0:
            v_r_ids = values2r_ids[(min_v, max_v)]
            if r_ids:
                v_r_ids = set(v_r_ids) & r_ids
            value = format_values(min_v, max_v)
            ess_count += len(v_r_ids)
            for r_id in sorted(v_r_ids):
                logging.info("%s\t%s: %s" % (value, r_id, model.reactions.get_by_id(r_id).build_reaction_string(True)))
            logging.info('')
    logging.info("==============================\nOTHER FLUXES\n==============================")
    for (min_v, max_v) in keys:
        if min_v * max_v == 0:
            v_r_ids = values2r_ids[(min_v, max_v)]
            if r_ids:
                v_r_ids = set(v_r_ids) & r_ids
            value = format_values(min_v, max_v)
            var_count += len(v_r_ids)
            for r_id in sorted(v_r_ids):
                logging.info("%s\t%s: %s" % (value, r_id, model.reactions.get_by_id(r_id).build_reaction_string(True)))
            logging.info('')
    logging.info("==============================")
    logging.info("In total, found %d essential reactions and %d variable reactions" % (ess_count, var_count))
    return r_id2bounds


def get_r_id2fva_bounds(model, threshold=None, rs=None):
    r_id2min_max = flux_variability_analysis(model, reaction_list=rs)
    r_id2bounds = {}
    for r_id, values in r_id2min_max.iteritems():
        min_v, max_v = round_value(values['minimum']), round_value(values['maximum'])
        if threshold is None or abs(min_v) > threshold or abs(max_v) > threshold:
            r_id2bounds[r_id] = min_v, max_v
    return r_id2bounds


def serialize_fva(model, r_id2bounds, path, r_ids=None, title=None):
    values2r_ids = defaultdict(set)
    for r_id, (min_v, max_v) in r_id2bounds.iteritems():
        values2r_ids[(min_v, max_v)].add(r_id)
    keys = sorted(values2r_ids.iterkeys(), key=lambda (min_v, max_v): (-abs(max_v - min_v), min_v, max_v))
    ess_count, var_count = 0, 0
    with open(path, 'w+') as f:
        if title:
            f.write(title + '\n')
        f.write("==============================\nESSENTIAL FLUXES\n==============================\n")
        for (min_v, max_v) in keys:
            if min_v * max_v > 0:
                v_r_ids = values2r_ids[(min_v, max_v)]
                if r_ids:
                    v_r_ids = set(v_r_ids) & r_ids
                value = format_values(min_v, max_v)
                ess_count += len(v_r_ids)
                for r_id in sorted(v_r_ids):
                    f.write("%s\t%s: %s\n" % (value, r_id,
                                              get_cobra_r_formula(model.reactions.get_by_id(r_id), comp=True)))
                f.write('\n')
        f.write("==============================\nOTHER FLUXES\n==============================\n")
        for (min_v, max_v) in keys:
            if min_v * max_v <= 0:
                v_r_ids = values2r_ids[(min_v, max_v)]
                if r_ids:
                    v_r_ids = set(v_r_ids) & r_ids
                value = format_values(min_v, max_v)
                var_count += len(v_r_ids)
                for r_id in sorted(v_r_ids):
                    f.write("%s\t%s: %s\n" % (value, r_id,
                                              get_cobra_r_formula(model.reactions.get_by_id(r_id), comp=True)))
                f.write('\n')
        f.write("==============================\n")
        f.write("In total, found %d essential reactions and %d variable reactions" % (ess_count, var_count))
    return ess_count, var_count


# The method find_blocked_reactions in COBRA fails as it tries to give unexpected parameters to the solver.
# This method is an analogue.
def get_blocked_reactions(cobra_model, the_reactions=None, tolerance_optimality=1e-9, open_exchanges=False, **kwargs):
    """Finds reactions that cannot carry a flux with the current
    exchange reaction settings for cobra_model, using flux variability
    analysis.
    """
    cobra_model = cobra_model.copy()
    if not the_reactions:
        the_reactions = cobra_model.reactions
    if open_exchanges:
        exchange_reactions = [x for x in cobra_model.reactions
                              if x.startswith('EX')]
        for the_reaction in exchange_reactions:
            if the_reaction.lower_bound >= 0:
                the_reaction.lower_bound = -1000
            if the_reaction.upper_bound >= 0:
                the_reaction.upper_bound = 1000
    flux_span_dict = flux_variability_analysis(cobra_model, the_reactions, **kwargs)
    return [k for k, v in flux_span_dict.items() if max(map(abs, v.values())) <= tolerance_optimality]


def constraint_reaction_of_interest(model, r_id, bound):
    reaction = model.reactions.get_by_id(r_id)
    up_b, l_b = abs(reaction.upper_bound), abs(reaction.lower_bound)
    set_bounds(reaction, bound, up_b, l_b)


def constraint_reactions_of_interest(model, rs, bound):
    for r_id in rs:
        constraint_reaction_of_interest(model, r_id, bound)


def format_r_name(r):
    return "%s: %s" % (r.id, r.name) if r.name.find(r.id) == -1 else r.name


def show_flux_variations_same_plot(model, figure_path, r_ids, r_id_var, bm_id, max_bound, log_x_scale=False,
                                   log_y_scale=False, remove_zeros=False, constrained=True, minimized=True):
    xs, r_id2ys = get_biomass_dependency(model, r_id_var, r_ids, bm_id, max_bound, remove_zeros, constrained, minimized)
    legend, data = [], []
    for r_id, ys in r_id2ys.iteritems():
        legend.append(format_r_name(model.reactions.get_by_id(r_id)))
        data.append(ys)
    r_var_name = format_r_name(model.reactions.get_by_id(r_id_var))
    create_plot(xs, data, legend, u'%s flux max bound (micromol/min/gDW)' % r_var_name,
                u'Flux values (micromol/min/gDW)', "", log_x_scale,
                log_y_scale)
    save_fig(figure_path)
    return xs, data, legend


def show_flux_variations_same_plot_separate_bm(model, figure_path, r_id_sets, r_id_var, bm_id, max_bound,
                                               log_x_scale=False,
                                               log_y_scale=False, remove_zeros=True,
                                               constrained=True, minimized=True):
    all_r_ids = {bm_id, r_id_var}
    for r_ids in r_id_sets:
        all_r_ids |= r_ids
    xs, r_id2ys = get_biomass_dependency(model, r_id_var, all_r_ids, bm_id, max_bound, remove_zeros, constrained,
                                         minimized)

    for r_id, ys in sorted(r_id2ys.iteritems(), key=lambda (r_id, ys): [abs(y) for y in ys]):
        logging.info("%s %s" % (ys, format_r_name(model.reactions.get_by_id(r_id))))

    r_var_name = format_r_name(model.reactions.get_by_id(r_id_var))
    num = len(r_id_sets)
    initialise_fig((16, 8 * num), 10)
    i = 1
    for r_ids in r_id_sets:
        legend, data = [], []
        for r_id in r_ids:
            ys = r_id2ys[r_id]
            legend.append(get_legend_part(model, r_id, ys))
            data.append(ys)
        create_subplot(xs, data, legend, "%s flux max bound" % r_var_name, "Fluxes", "Varying %s" % r_var_name,
                       log_x_scale,
                       log_y_scale, plot_num=200 + num * 10 + i)
        i += 2
    save_fig(figure_path)


def show_flux_variations_same_plot_separate_bms(model, figure_path, bm_id2r_ids, r_id_var, max_bound,
                                                log_x_scale=False,
                                                log_y_scale=False, remove_zeros=True,
                                                constrained=True, minimized=True):
    num = len(bm_id2r_ids)
    i = 1
    r_var_name = format_r_name(model.reactions.get_by_id(r_id_var))
    initialise_fig((16, 8 * num), 10)

    for bm_id, r_ids in bm_id2r_ids.iteritems():
        xs, r_id2ys = get_biomass_dependency(model, r_id_var, r_ids, bm_id, max_bound, remove_zeros, constrained,
                                             minimized)

        legend, data = [], []
        for r_id in r_ids:
            ys = r_id2ys[r_id]
            legend.append(get_legend_part(model, r_id, ys))
            data.append(ys)
        create_subplot(xs, data, legend, "%s flux max bound" % r_var_name, "Fluxes", "Varying %s" % r_var_name,
                       log_x_scale,
                       log_y_scale, plot_num=200 + num * 10 + i)
        i += 2
    save_fig(figure_path)


def minimize_fluxes(model, bm_id):
    biomass_value = model.solution.f
    r_id2min_max = flux_variability_analysis(model, allow_loops=False)
    r_ids = [r_id for (r_id, values) in r_id2min_max.iteritems() if values['minimum'] != values['maximum']]
    biomass_rn = model.reactions.get_by_id(bm_id)
    bm_l_b = biomass_rn.lower_bound
    bm_u_b = biomass_rn.upper_bound
    biomass_rn.lower_bound = biomass_value
    biomass_rn.upper_bound = biomass_value
    m2stoich = {}
    for r_id in r_ids:
        value = model.solution.x_dict[r_id]
        rn = model.reactions.get_by_id(r_id)
        for m in (rn.products if value > 0 else rn.reactants):
            if not m in m2stoich:
                m2stoich[m] = rn.get_coefficient(m)
            else:
                m2stoich[m] += rn.get_coefficient(m)
    fake_rn = Reaction("fake")
    fake_rn.add_metabolites(m2stoich)
    fake_rn.lower_bound = -1000
    fake_rn.upper_bound = 1000
    model.add_reaction(fake_rn)
    model.change_objective([model.reactions.get_by_id(bm_id)])
    model.optimize()
    model.change_objective([fake_rn])
    model.optimize(objective_sense='minimize')
    biomass_rn.lower_bound = bm_l_b
    biomass_rn.upper_bound = bm_u_b
    fake_rn.remove_from_model(model)


def get_legend_part(model, r_id, ys):
    r_name = format_r_name(model.reactions.get_by_id(r_id))
    m, M = min(ys), max(ys)
    value = format_values(m, M)
    return "%s %s" % (r_name, value)


def format_values(m, M):
    return "%g:\t%g" % (m, M) if m != M else "\tCONST\t%g" % m


def show_reaction_variations_same_plot(model, inhibited_rn_id, r_ids, bm_id, reaction_activity_rate, file_path, x_label,
                                       y_label, title, max_bound, log_x_scale=False, log_y_scale=False,
                                       remove_zeros=False,
                                       constrained=False, minimized=False, num_colors=0):
    if inhibited_rn_id:
        model.change_objective([model.reactions.get_by_id(bm_id)])
        model.optimize()
        bound = abs(model.solution.x_dict[inhibited_rn_id] * reaction_activity_rate)
        constraint_reaction_of_interest(model, inhibited_rn_id, bound)
    legend, data, xs = [], [], []
    for r_id in r_ids:
        xs, r_id2ys = get_biomass_dependency(model, r_id, {bm_id}, bm_id, max_bound, remove_zeros, constrained,
                                             minimized)
        ys = r_id2ys[bm_id]
        data.append(ys)
        r = model.reactions.get_by_id(r_id)
        legend.append(r.name.replace(r.id, '').strip())#format_r_name(r))
    initialise_fig((6, 6), 14)
    create_subplot(xs, data, legend, x_label, y_label, title, log_x_scale, log_y_scale, 111, legend_loc='upper left',
                   bb_to_anchor=(0, 1), num_colors=num_colors)
    save_fig(file_path)
    return xs, data, legend


def show_reaction_variations_same_plot_diff_biomass(model, r_ids, bm_ids, file_path, x_label,
                                                    y_label, title, max_bound, log_x_scale=False, log_y_scale=False,
                                                    remove_zeros=False,
                                                    constrained=False, minimized=False):
    legend, data, xs = [], [], []
    for r_id in r_ids:
        for bm_id in bm_ids:
            xs, r_id2ys = get_biomass_dependency(model, r_id, {bm_id}, bm_id, max_bound, remove_zeros, constrained,
                                                 minimized)
            ys = r_id2ys[bm_id]
            data.append(ys)
            legend.append(get_legend_part(model, bm_id, ys))
    create_plot(xs, data, legend, x_label, y_label, title, log_x_scale, log_y_scale)
    save_fig(file_path)


def show_reaction_values(model, var_r_id, r_ids, bm_id, file_path, x_label, y_label, title, max_bound,
                         log_x_scale=False, log_y_scale=False, remove_zeros=False, constrained=False, minimized=False):
    xs, r_id2ys = get_biomass_dependency(model, var_r_id, r_ids, bm_id, max_bound, remove_zeros, constrained,
                                         minimized)
    legend, data = [], []
    for r_id, ys in r_id2ys.iteritems():
        legend.append(format_r_name(model.reactions.get_by_id(r_id)))
        data.append(ys)
    initialise_fig((6, 6), 14)
    create_subplot(xs, data, legend, x_label, y_label, title, log_x_scale, log_y_scale, 111, legend_loc='upper left',
                   bb_to_anchor=(0, 1))
    save_fig(file_path)
    return xs, data, legend


def show_reaction_variations_same_plot_different_bounds(model, inhibited_rn_id, r_ids, bm_id, reaction_activity_rate,
                                                        file_path, x_label, y_label, title, max_bounds,
                                                        log_x_scale=False, log_y_scale=False, remove_zeros=False,
                                                        constrained=False, minimized=False):
    model.change_objective([model.reactions.get_by_id(bm_id)])
    model.optimize()
    bound = abs(model.solution.x_dict[inhibited_rn_id] * reaction_activity_rate)
    constraint_reaction_of_interest(model, inhibited_rn_id, bound)
    legend, data, x_data = [], [], []
    max_bounds = iter(max_bounds)
    for r_id in r_ids:
        max_bound = next(max_bounds)
        xs, r_id2ys = get_biomass_dependency(model, r_id, {bm_id}, bm_id, max_bound, remove_zeros, constrained,
                                             minimized)
        ys = r_id2ys[bm_id]
        data.append(ys)
        x_data.append(xs)
        legend.append(get_legend_part(model, r_id, ys))
    initialise_fig((6, 6), 10)
    create_subplot(x_data, data, legend, x_label, y_label, title, log_x_scale, log_y_scale, 111,
                   legend_loc='upper left', bb_to_anchor=(0, 1))
    save_fig(file_path)
