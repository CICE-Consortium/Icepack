## Overview
[![Build Status](https://travis-ci.org/CICE-Consortium/Icepack.svg?branch=master)](https://travis-ci.org/CICE-Consortium/Icepack)

This repository contains files describing the column physics of the sea ice model CICE, which is now maintained by the CICE Consortium.  For testing purposes and guidance for including Icepack in other sea ice host models, this repository also includes a driver and basic test suite.

Icepack is included in CICE as a git submodule.  Because Icepack is a submodule of CICE, Icepack and CICE development are handled independently with respect to the github repositories even though development and testing can be done together. 

## Obtaining Icepack

If you expect to make any changes to the code, we recommend that you first fork the Icepack repositories.  Basic instructions for working with CICE and Icepack are found in the Git Workflow Guidance, linked below and from the wiki      
https://github.com/CICE-Consortium/Icepack/wiki

Icepack may be obtained in several different ways:  [not yet tested]    
1.  clone the full repository    
See [Git Workflow Guidance](https://github.com/CICE-Consortium/About-Us/wiki/Git-Workflow-Guidance)    
2.  check out only a particular branch, version or tag    
In the workflow for step 1 above, substitute    
    git clone -b branch_name https://github.com/CICE-Consortium/Icepack.git local_directory_name   
or use svn    
   svn co https://github.com/CICE-Consortium/Icepack/branch_name    
where "branch name" can also be a version name    
3.  download a tarball for a particular version from the git releases: https://github.com/CICE-Consortium/Icepack/releases

## Documentation

Detailed and searchable online documentation of Icepack can be found at https://readthedocs.org/projects/cice-consortium-icepack/ . 

This site has the most up-to-date [HTML](http://cice-consortium-icepack.readthedocs.io/en/master/) and [PDF](https://media.readthedocs.org/pdf/cice-consortium-icepack/master/cice-consortium-icepack.pdf) living documentation from the master branch of the CICE-Consortium repository that will be updated regularly with code development. 

This site also has static documentation from each Icepack release.

More information about Icepack documentation can be found on the [Icepack Documentation Wiki page](https://github.com/CICE-Consortium/Icepack/wiki/Icepack-Documentation).

## More Information

The [Icepack wiki](https://github.com/CICE-Consortium/Icepack/wiki) page contains links to additional information, e.g.    
- complete documentation 
- larger files such as the gx1 grid, land mask, and forcing files
- testing data

The [Test-Results wiki](https://github.com/CICE-Consortium/Test-Results/wiki) has test results for both CICE and Icepack.

The [About-Us repository](https://github.com/CICE-Consortium/About-Us) includes background and supporting information about the CICE Consortium, including how to interact with it.   

See also our [FAQ](https://github.com/CICE-Consortium/About-Us/wiki/FAQ).   
