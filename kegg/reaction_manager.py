import logging
import urllib2
import datetime
import openpyxl
from xlsx_helper import save_data, HEADER_STYLE, add_values


def get_rns_by_elements(elements):
    rns = set()
    reactions_by_compounds = urllib2.urlopen(
        'http://rest.kegg.jp/link/reaction/%s' % "+".join(elements)).read()
    for line in reactions_by_compounds.split("\n"):
        if line.find("rn:") != -1:
            _, rn = line.split("\t")
            rns.add(rn)
    return rns


def get_compounds_by_rn(rn):
    compounds_by_reaction = urllib2.urlopen('http://rest.kegg.jp/link/compound/%s' % rn).read()
    return {c for c in {c.replace("%s\t" % rn, '').strip() for c in compounds_by_reaction.split("\n")} if
            c.find('cpd:') != -1}


def get_compounds_by_rp(rp):
    compounds_by_rp = urllib2.urlopen('http://rest.kegg.jp/link/compound/%s' % rp).read()
    return {c for c in {c.replace("%s\t" % rp, '').strip() for c in compounds_by_rp.split("\n")} if
            c.find('cpd:') != -1}


def get_rpairs_by_rn(rn):
    rpairs_by_reaction = urllib2.urlopen('http://rest.kegg.jp/link/rpair/%s' % rn).read()
    return {c for c in {c.replace("%s\t" % rn, '').strip() for c in rpairs_by_reaction.split("\n")} if
            c.find('rp:') != -1}


def get_rpairs_by_compounds(compounds):
    rpairs_by_compounds = urllib2.urlopen('http://rest.kegg.jp/link/rpair/%s' % '+'.join(compounds)).read()
    rpairs = set()
    for line in rpairs_by_compounds.split("\n"):
        if line.find("rp:") != -1:
            _, rp = line.split("\t")
            rpairs.add(rp)
    return rpairs


def get_pw_name(org, pw):
    result = urllib2.urlopen('http://rest.kegg.jp/find/pathway/%s' % p_id_specific2generic(org, pw)).read()
    result = result.replace("\n", '')
    pw, pw_name = result.split('\t')
    return pw_name


def get_rns():
    rn2formula = {}
    reactions = urllib2.urlopen('http://rest.kegg.jp/list/reaction').read()
    for line in reactions.split("\n"):
        if line.find("rn:") != -1:
            rn, r_formula = line.split("\t")
            rn2formula[rn] = r_formula
    return rn2formula


def get_kegg_r_info(rn):
    name, definition, equation, ec, orthology = '', '', '', '', ''
    try:
        reactions = urllib2.urlopen('http://rest.kegg.jp/get/%s' % rn).read()
        for line in reactions.split("\n"):
            if line.find("DEFINITION") != -1:
                definition = line.replace("DEFINITION", '').strip()
            elif line.find("NAME") != -1:
                name = line.replace("NAME", '').strip()
            elif line.find("EQUATION") != -1:
                equation = line.replace("EQUATION", '').strip()
            elif line.find("ENZYME") != -1:
                ec = line.replace("ENZYME", '').strip()
            elif line.find("ORTHOLOGY") != -1:
                orthology = line.replace("ORTHOLOGY", '').strip()
    except:
        pass
    return name, definition, equation, ec, orthology


def serialize_reactions(path):
    wb = openpyxl.Workbook()
    ws = wb.create_sheet(0, "Reactions")
    add_values(ws, 1, 1, ["Id", "Name", "Definition", "Formula", "EC", "Orthology"], HEADER_STYLE)
    row = 2
    rns = {it.split('\t')[0] for it in urllib2.urlopen('http://rest.kegg.jp/list/reaction').read().split("\n") if
             it.find('rn:') != -1}
    print (rns)
    for rn in rns:
        name, definition, equation, ec, orthology = get_kegg_r_info(rn)
        add_values(ws, row, 1, [rn.replace('rn:', '').strip(), name, definition, equation, ec, orthology])
        row += 1
    wb.save(path)


# serialize_reactions("/home/anna/Documents/IBGC/Models/KEGG_reactions.xlsx")

# print(get_kegg_r_info('R01513'))


def get_cpd2rps():
    cpd2rp = {}
    rps = urllib2.urlopen('http://rest.kegg.jp/list/rpair').read()
    for line in rps.split("\n"):
        if line.find("rp:") != -1:
            rp, info = line.split("\t")
            cpd2rp[info] = rp
    return cpd2rp


def get_relevant_rns(rpairs):
    rns = get_rns_by_elements(rpairs)
    for rn in rns:
        if rn:
            rn_rps = get_rpairs_by_rn(rn)
            if not rpairs - rn_rps:
                yield rn


def get_relevant_rpairs(compounds):
    rps = get_rpairs_by_compounds(compounds)
    for rp in rps:
        rp_cmps = get_compounds_by_rp(rp)
        if not rp_cmps - compounds:
            yield rp


def get_rns_by_compounds(r_compounds, p_compounds, cpd2rps=None):
    if not cpd2rps:
        cpd2rps = get_cpd2rps()
    comps = {"cpd:%s" % it for it in r_compounds} | {"cpd:%s" % it for it in p_compounds}
    if not r_compounds or not p_compounds:
        return iter([])
    rp = None
    for r in r_compounds:
        for p in p_compounds:
            cpd = "%s_%s" % (r, p) if r < p else (p, r)
            if cpd in cpd2rps:
                rp = cpd2rps[cpd]
                break
        if rp:
            break
    print(rp)
    return (rn for rn in get_relevant_rns({rp}) if get_compounds_by_rn(rn) == comps) if rp else iter([])
