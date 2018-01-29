##########################################################################
# @author: yann.leguen@cea.fr
# Copyright (C) CEA, 2017
# Distributed under the terms of the CeCILL-B license, as published by
# the CEA-CNRS-INRIA. Refer to the LICENSE file or to
# http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html for details.
##########################################################################

proc makeped {working_dir ped} {
    cd $working_dir
    exec mkdir -p pedigree
    cd pedigree
    field id id
    load pedigree $ped
    cd ..
}
