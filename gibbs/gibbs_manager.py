import os
import math
from decimal import Decimal
import re

import gibbs


DATA_REACTION_PH7 = "%s/data/kegg_reactions_CC_ph7.0.csv" % os.path.dirname(os.path.abspath(gibbs.__file__))
DATA_COMPOUND_PH7 = "%s/data/kegg_compounds_Alberty_ph7.0.csv" % os.path.dirname(os.path.abspath(gibbs.__file__))
R_PH2DATA = {7: DATA_REACTION_PH7}
C_PH2DATA = {7: DATA_COMPOUND_PH7}

# gas constant [kJ*mol^-1*K^-1]
R = Decimal(0.0083144627)

__author__ = 'anna'


def get_r_kegg2equilibrium_const(ph=7):
    if ph not in R_PH2DATA:
        raise ValueError('''
        Unfortunately, we don't have information for pH %d.
        You can try the following pH values instead: %s''' % (ph, ",".join([str(it) for it in R_PH2DATA.iterkeys()])))
    kegg2ke = {}
    with open(R_PH2DATA[ph], 'r') as f:
        # kegg, dG0_prime (kJ/mol), sigma[dG0] (kJ/mol), pH, I (mM), T (Kelvin), Note
        for line in f:
            # header
            if line.find("!") != -1:
                continue
            line = line.strip('\n').strip()
            kegg, g0, sigma, pH, I, T, note = line.split(',')
            if g0.strip():
                g0 = Decimal(g0.strip())
                kegg2ke[kegg.strip()] = float(get_ke(g0, Decimal(T.strip())))
    return kegg2ke


def get_c_kegg2g0(ph=7):
    if ph not in R_PH2DATA:
        raise ValueError('''
        Unfortunately, we don't have information for pH %d.
        You can try the following pH values instead: %s''' % (ph, ",".join([str(it) for it in C_PH2DATA.iterkeys()])))
    kegg2g0 = {}
    with open(C_PH2DATA[ph], 'r') as f:
        # kegg, Name, dG0_prime (kJ/mol), pH, I (mM), T (Kelvin), Note
        for line in f:
            line = line.strip('\n').strip()
            # header
            if line.find("!") != -1 or not line:
                continue
            # the name can contain commas (in that case it's in quotes).
            # as we are not interested by name, let's replace it by an empty string,
            # so it's possible commas do not mess with the delimiters:
            # 'C00198,"D-Glucono-1,5-lactone",-498.0,7.0,0.1,298.15,' --> 'C00198,,-498.0,7.0,0.1,298.15,'
            line = re.sub(r'\"(.+?)\"', "", line)
            kegg, name, g0, pH, I, T, note = line.split(',')
            if g0.strip():
                kegg2g0[kegg.strip()] = Decimal(g0.strip())
    return kegg2g0


def get_g0_equilibrium_const(kegg_compound2stoich, kegg_compound2g0):
    r_g0 = 0
    for (kegg, st) in kegg_compound2stoich.iteritems():
        if kegg not in kegg_compound2g0:
            return None, None
        r_g0 += kegg_compound2g0[kegg] * Decimal(st)
    return r_g0, get_ke(r_g0)


def get_ke(g0, t=Decimal(298.15)):
    return pow(Decimal(math.e), (-g0 / t) / R)
