FROM centos:8

# Set an encoding to make things work smoothly.
ENV LANG en_US.UTF-8
ENV PYTHONUNBUFFERED 1

RUN yum install -y mesa-libGLU-devel && \
    yum clean all && \
    useradd -m eman -G wheel && \
    echo '%wheel ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers && \
    echo 'Defaults:%wheel !requiretty' >> /etc/sudoers

USER eman
WORKDIR /home/eman

RUN curl -v -L https://cryoem.bcm.edu/cryoem/static/software/release-2.91/eman2.91_sphire1.4_sparx.linux64.sh -o eman.sh && \
    export EMAN_INSTALL_DONT_UPDATE_DEPS=1 && bash eman.sh -bp /home/eman/eman2_sphire_sparx && \
    rm eman.sh && \
    source /home/eman/eman2_sphire_sparx/etc/profile.d/conda.sh && \
    conda install eman-deps=26 -c cryoem -c defaults -c conda-forge && \
    conda clean --all --yes && \
    conda init

CMD [ "/bin/bash" ]

