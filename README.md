

#1 phenotype_extraction.py
Extract the median (or mean) parameter estimate (or z-stat) that
will be used as proxy of the activations in each areal

#2 heritability_analysis.py
run SOLAR heritability analysis of the phenotypes considered (tfMRI activations per areal)

#3 map_heritability.py
Map the heritability estimates and associated p-values on the template using Anatomist viewer

#4 heritability_cognitive_traits
run SOLAR heritability analysis on the cognitive scores provided by HCP and task fMRI performances

#5 bivariate_genetic_analysis.py
SOLAR bivariate genetic analysis of the pair of phenotypes (activation, cognitive score)

#6 map_bivariate_analysis.py
Map the phenotypic correlation and associated p-values on the template using Anatomist viewer
Map the shared genetic variance and associated p-values on the template using Anatomist viewer


makeped.tcl
SOLAR script to load HCP pedigree into SOLAR

pheno_analysis.tcl
SOLAR script to run the heritability analysis called from heritability_analysis.py

bivariate_analysis.tcl
SOLAR script to run the bivariate genetic analysis called from bivariate_genetic_analysis.py

equivalent_power_calculation.R
call the pedigree function from the R package RSP in order to get statistics from our dataset
NB: statistical power calculation in a pedigree is a difficult problem (cf Blangero et al 2013),
RSP (Raffa & Thompson 2016) provides a summary statistics allowing to compare the statistical power between sample pedigree.