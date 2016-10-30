
[![License](http://img.shields.io/:license-mit-blue.svg)](http://doge.mit-license.org)


# cmplx - Epidemic source detection library

## Compiling the library

First install the following dependencies. Those include dependencies to compile igraph locally, and opempi dependencies used by the library.

```bash
sudo apt-get install build-essential cmake bzr libtool libxml2-dev openmpi-common libopenmpi-dev
openmpi-bin
```

To download and compile igraph locally, run the script [```install_igraph.sh```](install_igraph.sh) in a new directory ``` igraph ```. It will download the latest version of igraph from Launchpad and compile it locally.

To enable linking, edit [```common/IgraphConfig.cmake```](common/IgraphConfig.cmake) and add local ```igraph/include``` dir path to ```IGRAPH_INCLUDE_DIR``` variable and local ```igraph/lib``` dir path to ```IGRAPH_LIBRARY``` variable.

Make sure the path to ```IgraphConfig.cmake``` is set properly in [```CMakeLists.txt```](CMakeLists.txt), [```common/CMakeLists.txt```](common/CMakeLists.txt) and [```simul/CMakeLists.txt```](simul/CMakeLists.txt).
Finally, create dir ```build``` and run 

```bash
cmake cmplx
make
```

###Compiling the tests
Tests are configured with ```cmake cmplx -Dtest=ON```.
Run them all with ```make test```.
Tests require [GTest](https://github.com/google/googletest) library. 

###Running
Run executables with ```mpiexec -n [#processes] executable [local flags]```.


