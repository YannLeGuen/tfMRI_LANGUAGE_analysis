##########################################################################
# @author: yann.leguen@cea.fr
# Copyright (C) CEA, 2017
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html for details.
##########################################################################

proc bivariate_analysis {work_dir out_dir trait trait2 phen_file} {

    field id IID
    outdir $out_dir
    cd $work_dir/pedigree
    load phenotypes $phen_file
    
    # Apply the inverse normal transformation
    define analyze_$trait = inormal_$trait 
    define analyze_$trait2 = inormal_$trait2
    # Define selected traits
    trait analyze_$trait analyze_$trait2
    
    # Add the covariates
    covariate age^1,2#sex
    covariate sex
    covariate group1

    # Previous options tested:
    # spormod 
    # house
    # -keephouse
    
    # The covariate educ might be removed because there might be shared
    # genetic variance effects, that we actually do not want to remove
    # when estimating heritability
    covariate educ

    # Perform the bivariate genetic analysis
    polygenic -testrhoe -testrhog -testrhop
}
