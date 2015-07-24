__author__ = 'anna'


URN_MIRIAM = "urn:miriam:"
IDENTIFIERS_ORG = "http://identifiers.org/"


def miriam_to_term_id(urn):
    urn = urn.strip()
    miriam_prefix = URN_MIRIAM
    identifiers_org_prefix = IDENTIFIERS_ORG
    if 0 == urn.find(miriam_prefix):
        urn = urn[len(miriam_prefix):]
        start = urn.find(":")
        if start != -1 and urn[start + 1:].find(":") != -1:
            urn = urn[start + 1:]
    elif 0 == urn.find(identifiers_org_prefix):
        urn = urn[len(identifiers_org_prefix):]
        start = urn.find("/")
        if start != -1:
            if urn[start + 1:].find(":") != -1:
                urn = urn[start + 1:]
            else:
                urn = urn.replace("/", ":")
    return urn.strip()


def to_identifiers_org_format(t_id, prefix=''):
    return "{0}{1}/{2}".format(IDENTIFIERS_ORG, prefix, t_id)