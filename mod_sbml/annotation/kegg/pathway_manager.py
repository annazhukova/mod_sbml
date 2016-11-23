from urllib.request import urlopen



ORG_HUMAN = "hsa"
ORG_SACE = "sce"
GENERIC = "map"


def p_id_specific2generic(org, pw):
    return pw.replace(org, "map")


def p_id_generic2specific(org, pw):
    return pw.replace("map", org)


def get_pw_by_reactions(rns, org):
    pws = set()
    rns = list(rns)
    i = 0
    while i < len(rns):
        j = min(i + 10, len(rns))
        pathways_by_reactions = urlopen(
            'http://rest.kegg.jp/link/pathway/%s+%s' % (org, "+".join(rns[i:j]))).read()
        for line in pathways_by_reactions.split("\n"):
            if line.find("path:map") != -1:
                rn, pw = line.split("\t")
                pws.add(pw)
        i += 10
    return pws


def get_reactions_by_pw(org, pw):
    reactions_by_pathway = urlopen('http://rest.kegg.jp/link/rn/%s+%s' % (org, pw)).read()
    return {r for r in {r.replace("%s\t" % pw, '').strip() for r in reactions_by_pathway.split("\n")} if
            r.find('rn:') != -1}


def get_pw_name(org, pw):
    try:
        result = urlopen('http://rest.kegg.jp/find/pathway/%s' % p_id_specific2generic(org, pw)).read()
    except:
        return None
    result = result.replace("\n", '')
    if result.find('\t') != -1:
        pw, pw_name = result.split('\t')
        return pw_name
    return None


def get_pw_by_organism(org):
    pw2name = {}
    pathways = urlopen('http://rest.kegg.jp/list/pathway/%s' % org).read()
    for line in pathways.split("\n"):
        if line.find("path:") != -1:
            pw, p_name = line.split("\t")
            pw2name[pw] = p_name
    return pw2name


def get_name2pw(org='map'):
    name2pw = {}
    pathways = urlopen('http://rest.kegg.jp/list/pathway/%s' % org).read().decode()
    for line in pathways.split("\n"):
        if line.find("path:") != -1:
            pw, p_name = line.split("\t")
            name2pw[p_name.lower().strip()] = pw[len('path:'):]
    return name2pw


def get_relevant_pathway_info(org, rns, threshold=0):
    pw2name = get_pw_by_organism(org)
    pws = get_pw_by_reactions(rns, org)

    pw2name_rs_ratio = {}

    for generic_pw in pws:
        pw = p_id_generic2specific(org, generic_pw)
        if pw not in pw2name:
            continue
        pw_rs = get_reactions_by_pw(org, generic_pw)
        if not pw_rs:
            continue
        rs = pw_rs & rns
        ratio = len(rs) * 1.0 / len(pw_rs)
        if ratio > threshold:
            pw2name_rs_ratio[pw] = pw2name[pw], rs, ratio

    return pw2name_rs_ratio
