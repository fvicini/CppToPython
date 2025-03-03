## Set Base Stage
FROM ubuntu:22.04 AS rbm-base

RUN apt update && DEBIAN_FRONTEND=noninteractive apt install -y \
apt-utils \
ssh \
software-properties-common \
apt-transport-https ca-certificates gnupg software-properties-common wget \
bash-completion

RUN add-apt-repository -y ppa:ubuntu-toolchain-r/test \
&& apt update \
&& DEBIAN_FRONTEND=noninteractive apt install -y \
gcc-13 g++-13 gcc-13-base gfortran-13 libopenmpi-dev  

RUN update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-13 100 \
&& update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-13 100 \
&& update-alternatives --install /usr/bin/gfortran gfortran /usr/bin/gfortran-13 100

RUN apt install -y python3.10 && \
apt install -y python3.10-distutils && \
apt install -y python3.10-dev && \
update-alternatives --install /usr/bin/python3 python3 /usr/bin/python3.10 100 && \
update-alternatives --install /usr/bin/python python /usr/bin/python3.10 100 && \
wget https://bootstrap.pypa.io/get-pip.py && \
python3 get-pip.py && \
rm -rf get-pip.py

RUN apt update && \
apt install -y cmake && \
apt install -y git && \
apt install -y cmake-curses-gui

RUN pip3 install jupyter --no-cache-dir && \
pip3 install numpy --no-cache-dir && \
pip3 install scipy --no-cache-dir && \
pip3 install matplotlib --no-cache-dir && \
pip3 install tqdm --no-cache-dir && \
pip3 install torch torchvision --index-url https://download.pytorch.org/whl/cpu --no-cache-dir

WORKDIR /content
RUN git clone https://github.com/fvicini/CppToPython.git

WORKDIR /content/CppToPython
RUN git submodule init
RUN git submodule update

WORKDIR /content/CppToPython/externals
RUN cmake -DINSTALL_VTK=OFF -DINSTALL_LAPACK=ON ../gedim/3rd_party_libraries && make -j4

WORKDIR /content/CppToPython/release
RUN cmake -DCMAKE_PREFIX_PATH="/content/CppToPython/externals/Main_Install/eigen3;/content/CppToPython/externals/Main_Install/triangle;/content/CppToPython/externals/Main_Install/tetgen;/content/CppToPython/externals/Main_Install/googletest;/content/CppToPython/externals/Main_Install/lapack" -DENABLE_TRIANGLE=ON -DENABLE_TETGEN=ON ../ && make -j4 GeDiM4Py

WORKDIR /content/
RUN git clone https://github.com/fvicini/RBM_ML_2025.git
WORKDIR /content/RBM_ML_2025

RUN echo "source /usr/share/bash-completion/completions/git;cd /content/RBM_ML_2025;git pull;" > /content/on_startup.sh
RUN chmod +x /content/on_startup.sh

CMD cd /content/CppToPython && git pull && \
cd /content/RBM_ML_2025 && git pull && \ 
jupyter notebook --ip 0.0.0.0 --port=8080 --no-browser --allow-root
