sudo apt-get install git
git clone https://github.com/imih/cmplx.git
sudo apt-get install build-essential cmake bzr libtool libxml2-dev openmpi-common libopenmpi-dev openmpi-bin
cd cmplx
./install_igraph.sh
cd ~/
git clone https://github.com/imih/master-thesis-data
ln -s master-thesis-data/network/ network
ln -s master-thesis-data/realizations realizations
mkdir ~/build
cd build
cmake ~/cmplx
make
cp ~/cmplx/pokreni.sh ~/build


