# cmplx - Epidemic source detection library

## Compiling the library

First install the following dependencies. Those include dependencies to compile igraph locally, and opempi dependencies used by the library.

```bash
sudo apt-get install build-essential cmake bzr libtool libxml2-dev openmpi-common libopenmpi-dev
openmpi-bin
```

To download and compile igraph locally, run the script [```install_igraph.sh```](install_igraph.sh) in a new directory ``` igraph ```. It will download the latest version of igraph from Launchpad and compile it locally.

To enable linking, edit [```common/IgraphConfig.cmake```](common/IgraphConfig.cmake) and add local ```igraph/include``` dir path to ```IGRAPH_INCLUDE_DIR``` variable and local ```igraph/lib``` dir path to ```IGRAPH_LIBRARY``` variable.

Make sure the path to ```IgraphConfig.cmake``` is set properly in [```CMakeLists.txt```](CMakeLists.txt), [```common/CMakeLists.txt```](common/CMakelists.txt) and [```CMakeLists.txt```](CMakeLists.txt).
Finally, create dir ```build``` and run 

```bash
cmake cmplx
make
```

###Compiling the tests



