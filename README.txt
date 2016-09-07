The scripts in this repository are NOT intended for re-use, will NOT be maintained
and are provided "as is" only. The intention is to illustrate the implementation 
of the rules-based and machine learning algorithms from Pesesky et al. 2016, and 
as a starting point for other groups to create their own implementations.

Contents:
The "GBASP_*.py" scripts are wrappers that can either executs the RB algorithm or
generate the input files for weka necessary for the LR algorithm. 
GBASP_resfams.py is heavily commented and intended for use with the 
E.coli_100.HMMAnnotation.txt example file. GBASP_resfinder.py and GBASP_card.py
have updated user interfaces, but are not commented, and are intended for use 
with the ResFinder_E.coli_100_matches.txt and CARD_E.coli_100_matches.txt
example annotation files, respectively.

The arff_converter.py script is used to convert the gene lists output by the 
GBASP_*.py scripts to arff format so they can be used as inputs to Weka 3 for
machine learning applications.

The E.coli_100_proteins.faa is an example protein sequence file for the example
isolate "E. coli 100" intended for use with all three scripts.

The Annotation.py and prediction_modules.py files are modules used by all three 
GBASP_*.py wrapper scripts.

The ReferenceFiles/ directory contains database specific keyword files used to 
classify resistance genes, as well as table files matching resistance profiles to
specific resistance gene variants.
