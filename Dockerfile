FROM docker.io/debian:bookworm-20220801
MAINTAINER Martin Heistermann <martin.heistermann@unibe.ch>

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
    binutils \
    ca-certificates \
    ccache \
    cmake \
    curl \
    g++-12 \
    gcc-12 \
    git \
    libboost-filesystem-dev \
    libboost-regex-dev \
    libboost-system-dev \
    libc-dev \
    libgmm++-dev\
    liblapack-dev \
    libopenblas64-serial-dev \
    libqt5gui5 \
    libtool \
    locales \
    make \
    ninja-build \
    qmake6 \
    qt5-qmake \
    qtbase5-dev \
    time \
    tzdata \
    wget

RUN sed -i '/en_US.UTF-8/s/^# //g' /etc/locale.gen && \
    dpkg-reconfigure --frontend=noninteractive locales && \
    update-locale LANG=en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US:en
ENV LC_ALL en_US.UTF-8

RUN ln -fs /usr/share/zoneinfo/Europe/Zurich /etc/localtime && dpkg-reconfigure -f noninteractive tzdata

RUN mkdir -p /opt && \
    cd /opt && \
    wget -q https://packages.gurobi.com/9.5/gurobi9.5.2_linux64.tar.gz && \
    echo 95d8ca18b7f86116ba834a27fd6228c5b1708ae67927e7ea0e954c09374a2d0f  gurobi9.5.2_linux64.tar.gz | sha256sum --check && \
    tar xf gurobi9.5.2_linux64.tar.gz && \
    ln -s gurobi952 gurobi

# gurobi's makefile calls g++ and is not easily overridable due to the special-character "C++"  variable name...
RUN ln -sf $(which g++-12) /usr/bin/g++
RUN ln -sf $(which gcc-12) /usr/bin/gcc

RUN make -C /opt/gurobi/linux64/src/build -j12 && \
    ln -sf /opt/gurobi/linux64/src/build/libgurobi_c++.a /opt/gurobi/linux64/lib/

COPY . /app

RUN mkdir /app/libs/CoMISo/build && \
          cd /app/libs/CoMISo/build && \
          cmake .. \
            -D CMAKE_C_COMPILER=gcc-12 \
            -D CMAKE_CXX_COMPILER=g++-12 \
            -D CMAKE_BUILD_TYPE=Release \
            -D CMAKE_CXX_FLAGS="-march=native" \
          && \
          cmake --build . -- -j12


RUN mkdir /app/build && \
    cd /app/build && \
    qmake ../quadwild && \
    make -j12 && \
    ln -s /app/build/quadwild /usr/local/bin/

RUN mkdir /app/components/quad_from_patches/build && \
    cd /app/components/quad_from_patches/build && \
    qmake .. && \
    make -j12 && \
    ln -s /app/components/quad_from_patches/build/quad_from_patches /usr/local/bin


ENV LD_LIBRARY_PATH="/app/libs/CoMISo/build/Build/lib/CoMISo:/opt/gurobi/linux64/lib/"


