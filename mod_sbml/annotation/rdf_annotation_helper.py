import libsbml

from mod_sbml.annotation.miriam_converter import miriam_to_term_id

__author__ = 'anna'


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
            if libsbml.BIOLOGICAL_QUALIFIER == term.getQualifierType() and qualifier_type == term.getBiologicalQualifierType():
                yield term


def get_qualifier_values(element, qualifier_type):
    for term in get_annotation_term_of_type(element, qualifier_type):
        for i in xrange(term.getNumResources()):
            yield term.getResourceURI(i).replace("%3A", ":")


def add_annotation(element, qualifier, annotation):
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
