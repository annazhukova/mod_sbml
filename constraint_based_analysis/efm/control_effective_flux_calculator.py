__author__ = 'anna'


def get_efm_efficiency(efm, r_id):
    r_id2coeff = efm.to_r_id2coeff()
    if not r_id2coeff or r_id not in r_id2coeff:
        return 0
    return r_id2coeff[r_id] / sum(abs(coeff) for coeff in r_id2coeff.itervalues())


def classify_efm_by_efficiency(id2efm, r_id):
    return {efm_id: get_efm_efficiency(efm, r_id) for (efm_id, efm) in id2efm.iteritems()}