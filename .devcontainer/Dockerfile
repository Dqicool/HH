FROM nvidia/cuda:11.1-devel-centos8
ARG USERNAME="qidong"
ARG GROUPNAME="zp"
ARG UID=1000
ARG GID=1000
WORKDIR /workdir
RUN yum -y update && yum -y install which git wget tar atlas-devel                                                      \
    redhat-lsb-core libX11-devel libXpm-devel libXft-devel libXext-devel                                                \
    openssl-devel glibc-devel rpm-build nano sudo python2 python38 python2-pip python38-pip &&                                                                    \
    wget https://cmake.org/files/v3.18/cmake-3.18.2-Linux-x86_64.tar.gz                                                 \
    && tar -C /usr/local --strip-components=1 --no-same-owner                                                           \
    -xvf cmake-*-Linux-x86_64.tar.gz     && rm cmake-*-Linux-x86_64.tar.gz
RUN echo '%wheel	ALL=(ALL)	NOPASSWD: ALL' >> /etc/sudoers &&                                                       \
    groupadd ${GROUPNAME} -g ${GID} &&                                                                                    \
    adduser ${USERNAME} -u ${UID} -G ${GROUPNAME} && chmod 755 /home/${USERNAME} &&                                                       \
    usermod -aG wheel ${USERNAME} &&                                                                                    \
    usermod -aG root ${USERNAME} &&                                                                                     \
    rm -rf /.dockerenv
ADD ./.globus   /home/${USERNAME}/.globus
ADD ./bashrc    /etc/bashrc
RUN chown -R ${UID}:${GID} /home/${USERNAME}/.globus && \
    chown 0:0 /etc/bashrc &&                            \
    ln -s /usr/bin/python2 /usr/bin/python &&           \
    ln -s /usr/bin/pip2 /usr/bin/pip 
USER ${USERNAME}
