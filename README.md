*This repository is under construction, and some resources listed here are not yet available*

## Overview

This repository contains files describing the column physics of the sea ice model CICE, which is now maintained by the CICE Consortium.  For testing purposes and guidance for including Icepack in other sea ice host models, this repository also includes a driver and basic test suite.

Icepack is included in CICE as a git submodule.  Because Icepack is a submodule of CICE, Icepack and CICE development are handled independently with respect to the github repositories even though development and testing can be done together. 

## Obtaining Icepack

If you expect to make any changes to the code, we recommend that you first fork the Icepack repositories.  Basic instructions for working with CICE and Icepack are found in the Git and Workflow Guide, linked from the wiki      
https://github.com/CICE-Consortium/Icepack/wiki

Icepack may be obtained in several different ways:  [not yet tested]    
1.  clone the full repository    
See [Git and Workflow Guide](https://docs.google.com/document/d/1rR6WAvZQT9iAMUp-m_HZ06AUCCI19mguFialsMCYs9o/edit?usp=sharing)    
2.  check out only a particular branch, version or tag    
In the workflow for step 1 above, substitute    
  [check this] git clone -b branch_name https://github.com/CICE-Consortium/Icepack.git local_directory_name   
or use svn    
  svn co https://github.com/CICE-Consortium/Icepack/branch_name    
where "branch name" can also be a version name    
3.  download a tarball for a particular version    
[how]

## More information

Detailed and searchable online documentation of Icepack can be found at https://cice-consortium.github.io/Icepack with the following sections:
- [Quick Start](https://cice-consortium.github.io/Icepack/intro/quickstart.html), which has instructions for running the model. 
- [Science guide](https://cice-consortium.github.io/Icepack/science_guide/index.html)
- [User guide](https://cice-consortium.github.io/Icepack/user_guide/index.html) with a [“Testing”](https://cice-consoritum.github.io/Icepack/user_guide/ug_testing.html) subsection with instructions for setting up standard tests (e.g. regression, restart).
- [Developer guide](https://cice-consortium.github.io/Icepack/developer_guide/index.html)
- [Index](https://cice-consortium.github.io/Icepack/icepack_index.html)
* Note that the documentation Table of Contents can be accessed at any point by clicking the upper left hand "Icepack 1.0 documentation" link.

The [Icepack wiki](https://github.com/CICE-Consortium/Icepack/wiki) page contains links to additional information, e.g.    
- complete documentation 
- larger files such as the gx1 grid, land mask, and forcing files
- testing data
- test results 

The ["About-Us" repository](https://github.com/CICE-Consortium/About-Us) includes background and supporting information about the CICE Consortium, including how to interact with it.    
