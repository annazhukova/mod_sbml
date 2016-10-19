import libsbml

__author__ = 'anna'

URN_MIRIAM = "urn:miriam:"
IDENTIFIERS_ORG = "http://identifiers.org/"


def miriam_to_term_id(urn):
    # for better web encoding colon is sometimes represented as %3A
    urn = normalise(urn)
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


def normalise(urn):
    return urn.strip().replace('%3A', ':').replace('%3a', ':')


def to_identifiers_org_format(t_id, prefix=''):
    if prefix and prefix[-1] != '/':
        prefix = '%s/' % prefix
    return "%s%s%s" % (IDENTIFIERS_ORG, prefix, t_id)


def get_is_qualifier():
    return libsbml.BQB_IS


def get_is_vo_qualifier():
    return libsbml.BQB_IS_VERSION_OF


def get_is_annotations(entity):
    return get_annotations(entity, get_is_qualifier())


def get_is_vo_annotations(entity):
    return get_annotations(entity, get_is_vo_qualifier())


def get_annotations(entity, qualifier):
    return (miriam_to_term_id(it) for it in get_qualifier_values(entity, qualifier))


def get_annotation_term_of_type(element, qualifier_type):
    cv_terms = element.getCVTerms()
    if cv_terms:
        for i in xrange(cv_terms.getSize()):
            term = cv_terms.get(i)
            if libsbml.BIOLOGICAL_QUALIFIER == term.getQualifierType() \
                    and qualifier_type == term.getBiologicalQualifierType():
                yield term


def get_qualifier_values(element, qualifier_type):
    for term in get_annotation_term_of_type(element, qualifier_type):
        for i in xrange(term.getNumResources()):
            yield normalise(term.getResourceURI(i))


def add_annotation(element, qualifier, annotation, prefix=''):
    annotation = to_identifiers_org_format(annotation, prefix)
    if not element.isSetMetaId():
        element.setMetaId("m_{0}".format(element.getId()))
    term = next(get_annotation_term_of_type(element, qualifier), None)
    if not term:
        term = libsbml.CVTerm(libsbml.BIOLOGICAL_QUALIFIER)
        term.setBiologicalQualifierType(qualifier)
        term.addResource(annotation)
        element.addCVTerm(term, True)
    else:
        term.addResource(annotation)


def remove_annotation(element, qualifier, annotation):
    if not element.isSetMetaId():
        element.setMetaId("m_{0}".format(element.getId()))
    for term in get_annotation_term_of_type(element, qualifier):
        term.removeResource(annotation)
