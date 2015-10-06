__author__ = 'anna'


def get_bounds(r, infinity=1000):
    k_law = r.getKineticLaw()
    r_lower_bound = -infinity if r.getReversible() else 0
    r_upper_bound = infinity
    if k_law:
        lb = k_law.getParameter("LOWER_BOUND")
        if lb:
            r_lower_bound = lb.getValue()
        ub = k_law.getParameter("UPPER_BOUND")
        if ub:
            r_upper_bound = ub.getValue()
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