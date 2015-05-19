from decimal import Decimal
import logging

from gibbs.gibbs_manager import get_r_kegg2equilibrium_const, get_c_kegg2g0, get_g0_equilibrium_const
from kegg.reaction_manager import get_rn2compounds
from sbml.sbml_manager import get_kegg_r_id, get_reactants, get_products, get_kegg_m_id


KEGG_PROTON = 'C00080'

__author__ = 'anna'

NO_KEGG_ID = 3
NO_KE = 2
NO_KEGG_INFO_FOUND = 1
DIFFERENT_FORMULAS = 2
SUCCESS = 0
REVERSED = 4


def compare_kegg_formulas(r, kegg2rs_ps, m2kegg):
    kegg = get_kegg_r_id(r)
    if not kegg:
        return NO_KEGG_ID

    if kegg not in kegg2rs_ps:
        return NO_KEGG_INFO_FOUND

    rs, ps = kegg2rs_ps[kegg]
    r_rs, r_ps = {m2kegg[m_id] if m_id in m2kegg else m_id for m_id in get_reactants(r)}, \
                 {m2kegg[m_id] if m_id in m2kegg else m_id for m_id in get_products(r)}
    rs -= {KEGG_PROTON}
    ps -= {KEGG_PROTON}
    r_rs -= {KEGG_PROTON}
    r_ps -= {KEGG_PROTON}

    if (r_rs != rs or r_ps != ps) and (r_rs != ps or r_ps != rs):
        return DIFFERENT_FORMULAS
    rev = len(ps & r_rs) + len(rs & r_ps) > len(rs & r_rs) + len(ps & r_ps)
    return REVERSED if rev else SUCCESS


def update_reaction_boundaries(model, m_id2kegg, infinity=1000, ph=7):
    """
    Update trivial reaction boundaries (0, +infinity or -infinity) in a given model based on equilibrium constant Ke,
    according to the following rule:
                    Ke  < 1.0e-9    =>   lower_bound = -1000, upper_bound = 0
        1.0e-9  <=  Ke  < 1.0e-6    =>   lower_bound = -1000, upper_bound = +10
        1.0e-6  <=  Ke  < 1.0e-3    =>   lower_bound = -1000, upper_bound = +100
        1.0e-3  <=  Ke  < 1.0e+3    =>   lower_bound = -1000, upper_bound = +1000
        1.0e+3  <=  Ke  < 1.0e+6    =>   lower_bound = -100,  upper_bound = +1000
        1.0e+6  <=  Ke  < 1.0e+9    =>   lower_bound = -10,   upper_bound = +1000
        1.0e+9  <=  Ke              =>   lower_bound = 0,     upper_bound = +1000.
    Non-trivial reaction boundaries (other than 0, +infinity or -infinity) are considered as manually curated
    and are not updated.

    :param model: cobra model to be updated
    :param m_id2kegg: dict, metabolite id to KEGG compound database (http://www.kegg.jp/kegg/compound) id,
    e.g. {"m_239": "C00430", "m_566": "C00931"}
    :param infinity: int, the value of a reaction bound that is considered as a positive infinity in the model,
    by default 1000
    :param ph: int, pH value for choosing the equilibrium constant, by default 7.
    """
    # kegg2ke = get_r_kegg2equilibrium_const(ph)
    # kegg2rs_ps = get_rn2compounds()
    # rs_ps2kegg = get_compounds2rn()

    m_id2kegg = {m.id: get_kegg_m_id(m) for m in model.getListOfSpecies()}

    kegg_compound2g0 = get_c_kegg2g0(ph)

    def get_ke(r):
        m_kegg2st = {}
        for m_id, st in get_reactants(r, True):
            kegg = m_id2kegg[m_id]
            if not kegg:
                return None, None
            if kegg in m_kegg2st:
                m_kegg2st[kegg] += -st
            else:
                m_kegg2st[kegg] = -st
        for m_id, st in get_products(r, True):
            kegg = m_id2kegg[m_id]
            if not kegg:
                return None, None
            if kegg in m_kegg2st:
                m_kegg2st[kegg] += st
            else:
                m_kegg2st[kegg] = st
        return get_g0_equilibrium_const(m_kegg2st, kegg_compound2g0)

    r_id2ke = {}
    for r in model.reactions:
        g0, ke = get_ke(r)
        update_boundaries_sbml(r, ke, infinity)
        r_id2ke[r.id] = g0, ke
    return r_id2ke


def update_boundaries(r, kegg2ke, r_id2kegg, m2kegg, kegg2rs_ps, infinity=1000):
    if not infinity:
        infinity = 1000
    infinity = abs(infinity)
    if r.id not in r_id2kegg:
        return NO_KEGG_ID, None
    kegg = r_id2kegg[r.id]

    # _, _, equation, _, _ = get_kegg_r_info('rn:%s' % kegg)
    # if not equation:
    #     return NO_KEGG_INFO_FOUND, None
    #
    # rs, ps = get_rs_ps_by_kegg_equation(equation)

    if kegg not in kegg2rs_ps:
        return NO_KEGG_INFO_FOUND, None

    rs, ps = kegg2rs_ps[kegg]
    r_rs, r_ps = {m2kegg[m.id] for m in r.reactants if m.id in m2kegg}, {m2kegg[m.id] for m in r.products if
                                                                         m.id in m2kegg}

    reversed = len(ps & r_rs) + len(rs & r_ps) > len(rs & r_rs) + len(ps & r_ps)

    if (r_rs != rs or r_ps != ps) and (r_rs != ps or r_ps != rs):
        return DIFFERENT_FORMULAS, reversed

    if reversed:
        logging.info(
            "Reaction %s is defined in the opposite direction to the corresponding KEGG reaction %s" % (r.id, kegg))

    if kegg not in kegg2ke:
        return NO_KE, None

    ke = kegg2ke[kegg]

    if ke < Decimal(1.0e3):
        # if the reaction is allowed in that direction,
        # let's keep their bound,
        # otherwise we put -infinity
        if not reversed:
            if r.lower_bound >= 0:
                r.lower_bound = -1000
        else:
            if r.upper_bound <= 0:
                r.upper_bound = 1000

        # if the reaction is allowed in that direction,
        # and the bound is not infinity (i.e. was manually specified),
        # let's keep their bound,
        # otherwise we update it according to Ke
        if not reversed and (r.upper_bound <= 0 or r.upper_bound >= infinity) \
                or reversed and (r.lower_bound >= 0 or r.lower_bound <= -infinity):
            if ke < Decimal(1.0e-9):
                if not reversed:
                    r.upper_bound = 0
                else:
                    r.lower_bound = 0
            elif Decimal(1.0e-9) <= ke < Decimal(1.0e-5):
                if not reversed:
                    r.upper_bound = 10
                else:
                    r.lower_bound = -10
            elif Decimal(1.0e-5) <= ke < Decimal(1.0e-3):
                if not reversed:
                    r.upper_bound = 100
                else:
                    r.lower_bound = -100
            else:
                if not reversed:
                    r.upper_bound = 1000
                else:
                    r.lower_bound = -1000
    else:
        # if the reaction is allowed in that direction,
        # let's keep their bound,
        # otherwise we put infinity
        if not reversed:
            if r.upper_bound <= 0:
                r.upper_bound = 1000
        else:
            if r.lower_bound >= 0:
                r.lower_bound = -1000

        # if the reaction is allowed in that direction,
        # and the bound is not -infinity (i.e. was manually specified),
        # let's keep their bound,
        # otherwise we update it according to Ke
        if not reversed and (r.lower_bound >= 0 or r.lower_bound <= -infinity) \
                or reversed and (r.upper_bound <= 0 or r.upper_bound >= infinity):
            if ke >= Decimal(1.0e9):
                if not reversed:
                    r.lower_bound = 0
                else:
                    r.upper_bound = 0
            elif Decimal(1.0e9) > ke >= Decimal(1.0e5):
                if not reversed:
                    r.lower_bound = -10
                else:
                    r.upper_bound = 10
            elif Decimal(1.0e5) > ke >= Decimal(1.0e3):
                if not reversed:
                    r.lower_bound = -100
                else:
                    r.upper_bound = 100
            else:
                if not reversed:
                    r.lower_bound = -1000
                else:
                    r.upper_bound = 1000
    return SUCCESS, reversed


def get_bounds(r, infinity=1000):
    k_law = r.getKineticLaw()
    if not k_law:
        k_law = r.createKineticLaw()
        k_law.setFormula("FLUX_VALUE")

    lb = k_law.getParameter("LOWER_BOUND")
    if not lb:
        lb = k_law.createParameter()
        lb.setId("LOWER_BOUND")
        lb.setValue(-infinity if r.getReversible() else 0)
        lb.setUnits("mumol_per_gDW_per_min")
    r_lower_bound = lb.getValue()

    ub = k_law.getParameter("UPPER_BOUND")
    if not ub:
        ub = k_law.createParameter()
        ub.setId("UPPER_BOUND")
        ub.setValue(infinity)
        ub.setUnits("mumol_per_gDW_per_min")
    r_upper_bound = ub.getValue()

    if not k_law.getParameter("FLUX_VALUE"):
        fv = k_law.createParameter()
        fv.setId("FLUX_VALUE")
        fv.setValue(0)
        fv.setUnits("mumol_per_gDW_per_min")

    if not k_law.getParameter("OBJECTIVE_COEFFICIENT"):
        oc = k_law.createParameter()
        oc.setId("OBJECTIVE_COEFFICIENT")
        ub.setValue(0)

    return r_lower_bound, r_upper_bound


def set_bounds(r, r_lower_bound, r_upper_bound):
    k_law = r.getKineticLaw()
    if not k_law:
        k_law = r.createKineticLaw()
        k_law.setFormula("FLUX_VALUE")

    if not k_law.getParameter("FLUX_VALUE"):
        fv = k_law.createParameter()
        fv.setId("FLUX_VALUE")
        fv.setValue(0)
        fv.setUnits("mumol_per_gDW_per_min")

    if not k_law.getParameter("OBJECTIVE_COEFFICIENT"):
        oc = k_law.createParameter()
        oc.setId("OBJECTIVE_COEFFICIENT")
        oc.setValue(0)

    lb = k_law.getParameter("LOWER_BOUND")
    if not lb:
        lb = k_law.createParameter()
        lb.setId("LOWER_BOUND")
        lb.setUnits("mumol_per_gDW_per_min")
    lb.setValue(r_lower_bound)

    ub = k_law.getParameter("UPPER_BOUND")
    if not ub:
        ub = k_law.createParameter()
        ub.setId("UPPER_BOUND")
        ub.setUnits("mumol_per_gDW_per_min")
    ub.setValue(r_upper_bound)


def get_better_r_bounds(r, ke, rev, infinity=1000):
    if not infinity:
        infinity = 1000
    infinity = abs(infinity)

    r_lower_bound, r_upper_bound = get_bounds(r, infinity)
    if ke < Decimal(1.0e3):
        # if the reaction is allowed in that direction,
        # let's keep their bound,
        # otherwise we put -infinity
        if not rev:
            if r_lower_bound >= 0:
                r_lower_bound = -infinity
        else:
            if r_upper_bound <= 0:
                r_upper_bound = infinity

        # if the reaction is allowed in that direction,
        # and the bound is not infinity (i.e. was manually specified),
        # let's keep their bound,
        # otherwise we update it according to Ke
        if not rev and (r_upper_bound <= 0 or r_upper_bound >= infinity) \
                or rev and (r_lower_bound >= 0 or r_lower_bound <= -infinity):
            if ke < Decimal(1.0e-9):
                if not rev:
                    r_upper_bound = 0
                else:
                    r_lower_bound = 0
            elif Decimal(1.0e-9) <= ke < Decimal(1.0e-6):
                if not rev:
                    r_upper_bound = infinity / 100
                else:
                    r_lower_bound = -infinity / 100
            elif Decimal(1.0e-5) <= ke < Decimal(1.0e-3):
                if not rev:
                    r_upper_bound = infinity / 10
                else:
                    r_lower_bound = -infinity / 10
            else:
                if not rev:
                    r_upper_bound = infinity
                else:
                    r_lower_bound = -infinity
    else:
        # if the reaction is allowed in that direction,
        # let's keep their bound,
        # otherwise we put infinity
        if not rev:
            if r_upper_bound <= 0:
                r_upper_bound = infinity
        else:
            if r_lower_bound >= 0:
                r_lower_bound = -infinity

        # if the reaction is allowed in that direction,
        # and the bound is not -infinity (i.e. was manually specified),
        # let's keep their bound,
        # otherwise we update it according to Ke
        if not rev and (r_lower_bound >= 0 or r_lower_bound <= -infinity) \
                or rev and (r_upper_bound <= 0 or r_upper_bound >= infinity):
            if ke >= Decimal(1.0e9):
                if not rev:
                    r_lower_bound = 0
                else:
                    r_upper_bound = 0
            elif Decimal(1.0e9) > ke >= Decimal(1.0e6):
                if not rev:
                    r_lower_bound = -infinity / 100
                else:
                    r_upper_bound = infinity / 100
            elif Decimal(1.0e5) > ke >= Decimal(1.0e3):
                if not rev:
                    r_lower_bound = -infinity / 10
                else:
                    r_upper_bound = infinity / 10
            else:
                if not rev:
                    r_lower_bound = -infinity
                else:
                    r_upper_bound = infinity
    return r_lower_bound, r_upper_bound


def update_boundaries_sbml(r, ke, rev, infinity=1000):
    if ke is None:
        return NO_KE
    r_lower_bound, r_upper_bound = get_better_r_bounds(r, ke, rev, infinity)
    set_bounds(r, r_lower_bound, r_upper_bound)
    return SUCCESS


def update_bounds(model, infinity=1000):
    kegg2ke = get_r_kegg2equilibrium_const(ph=7)
    kegg2rs_ps = get_rn2compounds()
    m2kegg = {m.id: get_kegg_m_id(m) for m in model.getListOfSpecies()}
    for r in model.getListOfReactions():
        kegg = get_kegg_r_id(r)
        if not kegg or kegg not in kegg2ke:
            continue
        result = compare_kegg_formulas(r, kegg2rs_ps, m2kegg)
        if result in [SUCCESS, REVERSED]:
            ke = kegg2ke[kegg]
            r_lower_bound, r_upper_bound = get_better_r_bounds(r, ke, result == REVERSED, infinity)
            set_bounds(r, r_lower_bound, r_upper_bound)