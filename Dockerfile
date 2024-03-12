FROM spack/ubuntu-jammy:0.21.2 as builder

RUN apt-get update -y && apt-get upgrade -y && apt-get install -y wget && \
apt-get clean -y && \
apt-get autoremove -y && \
apt-get autoclean -y && \
rm -rf /var/lib/apt/lists/*

RUN wget --progress=bar:force:noscroll https://gitlab.lrz.de/hpcsoftware/Peano/-/archive/p4/Peano-p4.tar.gz \
&& tar -xvf Peano-p4.tar.gz \
&& spack repo add Peano-p4/spack

RUN mkdir /opt/spack-environment \
&&  (echo spack: \
&&   echo '  packages:' \
&&   echo '    all:' \
&&   echo '      target: [x86_64]' \
&&   echo '' \
&&   echo '  specs:' \
&&   echo '  - exahype2+omp' \
&&   echo '  - py-pip' \
&&   echo '  - py-jupyterlab' \
&&   echo '  - paraview+ipo+python ^llvm ~clang ~flang ~lldb ~lld libcxx=none ~gold +llvm_dylib' \
&&   echo '' \
&&   echo '  concretizer:' \
&&   echo '    unify: true' \
&&   echo '  config:' \
&&   echo '    install_tree: /opt/software' \
&&   echo '  view: /opt/views/view') > /opt/spack-environment/spack.yaml

RUN cd /opt/spack-environment && spack env activate . && spack install --fail-fast && spack gc -y
RUN cd /opt/spack-environment && spack env activate --sh -d . > activate.sh

FROM peanoframework/base:latest
LABEL peano-framework.org=http://www.peano-framework.org/

ENV FORCE_UNSAFE_CONFIGURE=1 \
DEBIAN_FRONTEND=noninteractive \
TZ=UTC \
DISPLAY=DISPLAY=:99.0 \
PYVISTA_OFF_SCREEN=true \
PYVISTA_USE_IPYVTK=true \
OMP_PLACES="cores" \
OMP_PROC_BIND="spread" \
OMPI_ALLOW_RUN_AS_ROOT=1 \
OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

RUN apt-get update -y && apt-get upgrade -y && apt-get install -y --no-install-recommends \
tini \
xvfb && \
apt-get clean -y && \
apt-get autoremove -y && \
apt-get autoclean -y && \
rm -rf /var/lib/apt/lists/*

RUN python3 -m pip install --upgrade \
pip \
panel \
pyvista \
ipywidgets \
vtk

COPY --from=builder /opt/spack-environment /opt/spack-environment
COPY --from=builder /opt/software /opt/software

# paths.view is a symlink, so copy the parent to avoid dereferencing and duplicating it
COPY --from=builder /opt/views /opt/views

RUN { \
      echo '#!/bin/sh' \
      && echo '.' /opt/spack-environment/activate.sh \
      && echo 'which Xvfb' \
      && echo 'Xvfb :99 -screen 0 1024x768x24 > /dev/null 2>&1 &' \
      && echo 'sleep 3' \
      && echo 'set +x' \
      && echo 'jupyter lab --allow-root --port=9999 --no-browser --ip=0.0.0.0' \
      && echo 'exec "$@"'; \
    } > /entrypoint.sh \
&& chmod a+x /entrypoint.sh \
&& ln -s /opt/views/view /opt/view

VOLUME ["/training"]
WORKDIR /training
ENTRYPOINT ["tini", "-s", "-g", "/entrypoint.sh", "--"]
