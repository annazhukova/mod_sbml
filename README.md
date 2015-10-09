# mod_sbml

mod_sbml implements utilities for working with metabolic models in SBML format [http://sbml.org], 
that include:
 - model annotation with identifiers from:
    * the KEGG database [http://www.genome.jp/kegg], 
    * the ChEBI Ontology [http://www.ebi.ac.uk/chebi/],
    * the 'cellular component' branch of the Gene Ontology [http://geneontology.org/];
 - model serialization to csv format [https://en.wikipedia.org/wiki/Comma-separated_values];
 - submodel extraction.
 
mod_sbml also includes utilities for parsing ontologies in OBO flat file format [https://github.com/owlcollab/oboformat].
 
mod_sbml uses the libsbml package [http://sbml.org/Software/libSBML] for manipulation of SBML files.

