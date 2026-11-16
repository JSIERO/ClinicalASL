# syntax=docker/dockerfile:1
# =========================
# ClinicalASL: MRI DIAMOX ASL CBF/AAT/ATA/CVR DICOM-in / DICOM-out container for IMAGR.
# Lean multi-stage build on Ubuntu 22.04, CPU only, runs offline (--network none).
#   dcm2niix : dcm2niix v1.0.20220720 with JPEG 2000 (OpenJPEG) + JPEG-LS (CharLS) for PACS single-frame DICOM
#   pyenv    : micromamba env "qi" = QASL (licensed, bitbucket) + CPU-only torch + HD-BET + ClinicalASL deps, cleaned
#   final    : ubuntu 22.04 runtime libs + the above, no compilers, no git, no micromamba, no QASL key
#
# Build (QASL key passed as BuildKit secret, never stored in a layer):
#   docker build --secret id=qasl_key,src=G:\DATA\MOYAMOYA\IMAGER\qi-key-jeroen \
#       --build-arg APP_VERSION=1.0.4 -t jeroensnl/mri_diamox_umcu_clinicalasl_cvr_imager:1.0.4 .
# =========================

# =========================
# dcm2niix (build stage)
# =========================
FROM ubuntu:22.04 AS dcm2niix
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates git cmake make g++ pkg-config zlib1g-dev libopenjp2-7-dev libcharls-dev \
    && apt-get clean && rm -rf /var/lib/apt/lists/*
RUN git clone --depth 1 --branch v1.0.20220720 https://github.com/rordenlab/dcm2niix.git /src \
    && mkdir /src/build && cd /src/build \
    && cmake -DCMAKE_BUILD_TYPE=Release -DUSE_OPENJPEG=ON -DUSE_CHARLS=ON .. \
    && make -j"$(nproc)" && strip bin/dcm2niix && cp bin/dcm2niix /usr/local/bin/ \
    && dcm2niix -v | head -n 3

# =========================
# Python environment (build stage)
# Same env as qasl_setup builds (micromamba env "qi", python 3.11, pip installs from QI bitbucket), but
#  - reproducible: QI packages at the exact commits of the validated production image (docker/qasl-commits.txt),
#    all PyPI packages constrained to the production versions (docker/qi-constraints.txt). qasl_setup itself only
#    installs HEAD (QASL 1.6.0rc3 at the time of writing, which needs the ANTs CLI; there is no release/1.4.2 branch).
#  - torch installed CPU-only first, so nothing pulls the CUDA wheels (nvidia-*, triton: ~3.4 GB).
# =========================
FROM ubuntu:22.04 AS pyenv
ENV DEBIAN_FRONTEND=noninteractive \
    MAMBA_ROOT_PREFIX=/opt/conda \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1
RUN apt-get update && apt-get install -y --no-install-recommends \
        ca-certificates git openssh-client wget unzip sed gcc g++ \
    && apt-get clean && rm -rf /var/lib/apt/lists/*
COPY --from=mambaorg/micromamba:1.5.10 /bin/micromamba /usr/local/bin/micromamba

ARG TORCH_VERSION=2.13.0
ARG TORCHVISION_VERSION=0.28.0
ARG QASL_VERSION=1.4.2
ARG Q=/opt/conda/envs/qi
COPY docker/qi-constraints.txt docker/qasl-commits.txt /tmp/
ENV PIP_CONSTRAINT=/tmp/qi-constraints.txt
RUN micromamba create -y -n qi -c conda-forge python=3.11 pip \
    && ${Q}/bin/pip install --index-url https://download.pytorch.org/whl/cpu --extra-index-url https://pypi.org/simple \
        "torch==${TORCH_VERSION}" "torchvision==${TORCHVISION_VERSION}" \
    && ${Q}/bin/pip install -r /tmp/qi-constraints.txt

# QASL packages, --no-deps: their deps are the frozen set above (production has numpy 2.0.1 although qasl_qc
# declares numpy<2). License key via BuildKit secret: never in a layer.
RUN --mount=type=secret,id=qasl_key,required=true \
    export GIT_SSH_COMMAND="ssh -i /run/secrets/qasl_key -o IdentityAgent=none -o StrictHostKeyChecking=no -o UserKnownHostsFile=/dev/null" \
    && grep -vE "^\s*(#|$)" /tmp/qasl-commits.txt | while read -r pkg commit; do \
        echo "== ${pkg} @ ${commit}" \
        && ${Q}/bin/pip install --no-deps "git+ssh://git@ssh.bitbucket.org/quantified-imaging/${pkg}@${commit}" || exit 1; \
    done \
    && ${Q}/bin/pip list 2>/dev/null | grep -E "^(qasl|ssvb|vaby|toblerone|regtricks) " \
    && ${Q}/bin/python -c "import importlib.metadata as m, sys; v=m.version('qasl'); print('qasl', v); sys.exit(0 if v.startswith('${QASL_VERSION}') else 1)"

# HD-BET v2.0.1 + model parameters (non-editable install, parameter folder patched to /opt/hd-bet_params)
RUN git clone --depth 1 --branch v2.0.1 https://github.com/MIC-DKFZ/HD-BET.git /tmp/HD-BET \
    && sed -i "s|folder_with_parameter_files = .*|folder_with_parameter_files = '/opt/hd-bet_params'|" \
        /tmp/HD-BET/HD_BET/paths.py \
    && ${Q}/bin/pip install --no-deps /tmp/HD-BET \
    && mkdir -p /opt/hd-bet_params \
    && wget -q https://zenodo.org/records/14445620/files/release_v1.5.0.zip -O /tmp/hd-bet.zip \
    && unzip -q /tmp/hd-bet.zip -d /opt/hd-bet_params \
    && rm -rf /tmp/hd-bet.zip /tmp/HD-BET

# ClinicalASL python deps (antspyx==0.5.4 pinned: 0.6.x registration performs badly)
COPY requirements.txt /tmp/requirements.txt
RUN ${Q}/bin/pip install -r /tmp/requirements.txt && (${Q}/bin/pip check || true)

# toblerone (used by ssvb) AOT-compiles a numba module on first import and needs a C compiler for that:
# do it here so the final image ships the compiled .so and needs no compiler
RUN cd ${Q}/lib/python3.11/site-packages/toblerone && ${Q}/bin/python numba_precomp_src.py \
    && ls numba_precomp*.so && ${Q}/bin/python -c "import toblerone.numba_precomp"

# Clean: dev/test packages from qasl_setup, caches, headers, static libs, test folders
RUN ${Q}/bin/pip uninstall -y pytest pytest-cov pytest-xdist flake8 cython qasl_system_tests 2>/dev/null; \
    if ${Q}/bin/pip list 2>/dev/null | grep -iE "^(nvidia-|triton |cuda-)"; then echo "CUDA wheels present" && exit 1; fi \
    && micromamba clean -a -y \
    && SP=${Q}/lib/python3.11/site-packages \
    && rm -rf ${Q}/include ${Q}/share/doc ${Q}/share/man ${Q}/conda-meta /opt/conda/pkgs /root/.cache \
        ${SP}/torch/include ${SP}/torch/share ${SP}/tensorflow/include \
    && find ${Q} -name "*.a" -delete \
    && find ${SP} -depth -type d \( -name tests -o -name test \) \
        -not -path "*/torch/*" -not -path "*/tensorflow/*" -not -path "*/numpy/*" -not -path "*/keras/*" \
        -exec rm -rf {} +

# =========================
# Final image
# =========================
FROM ubuntu:22.04
ARG APP_VERSION=1.0.4
LABEL org.opencontainers.image.title="ClinicalASL" \
      org.opencontainers.image.version="${APP_VERSION}" \
      org.opencontainers.image.authors="Jeroen Siero, UMC Utrecht" \
      org.opencontainers.image.source="https://github.com/JSIERO/ClinicalASL" \
      org.opencontainers.image.description="ASL CBF/AAT/ATA/CVR (QASL) for IMAGR - not CE marked / FDA cleared"

ENV DEBIAN_FRONTEND=noninteractive \
    LANG=C.UTF-8 LC_ALL=C.UTF-8 \
    PATH=/opt/conda/envs/qi/bin:/usr/local/sbin:/usr/local/bin:/usr/sbin:/usr/bin:/sbin:/bin \
    MPLBACKEND=Agg \
    MKL_THREADING_LAYER=GNU \
    TF_CPP_MIN_LOG_LEVEL=2

# =========================
# System runtime libraries only (no compilers)
# =========================
RUN apt-get update && apt-get install -y --no-install-recommends \
        libgomp1 zlib1g libopenjp2-7 libcharls2 \
    && apt-get clean && rm -rf /var/lib/apt/lists/*

# =========================
# Python env, HD-BET parameters, dcm2niix (same paths as the previous image)
# =========================
COPY --from=pyenv /opt/conda/envs/qi /opt/conda/envs/qi
COPY --from=pyenv /opt/hd-bet_params /opt/hd-bet_params
COPY --from=dcm2niix /usr/local/bin/dcm2niix /opt/conda/envs/qi/bin/dcm2niix

# =========================
# ClinicalASL pipeline + entrypoint
# =========================
COPY python /opt/ClinicalASL/python
COPY docker-entrypoint.sh /usr/local/bin/docker-entrypoint.sh
RUN sed -i "s/^__version__ = .*/__version__ = \"${APP_VERSION}\"/" /opt/ClinicalASL/python/clinical_asl_pipeline/__version__.py \
    && sed -i 's/\r$//' /usr/local/bin/docker-entrypoint.sh && chmod +x /usr/local/bin/docker-entrypoint.sh

# =========================
# FINAL CHECK (fails the build if anything is missing)
# =========================
RUN python -c "import numpy, scipy, nibabel, pydicom, matplotlib, ants, torch, HD_BET, tensorflow, qasl, ssvb, vaby; \
assert torch.version.cuda is None, 'torch is not the CPU build'; \
print('FINAL OK numpy', numpy.__version__, 'ants', ants.__version__, 'torch', torch.__version__, \
'tf', tensorflow.__version__, 'pydicom', pydicom.__version__)" \
    && ! ldd /opt/conda/envs/qi/bin/dcm2niix | grep -q "not found" && dcm2niix -v | head -n 1 \
    && qasl --help >/dev/null && hd-bet -h >/dev/null && echo "QASL + HD-BET OK" \
    && test -d /opt/hd-bet_params \
    && cd /opt/ClinicalASL/python \
    && python -c "from clinical_asl_pipeline import main_pipeline; from clinical_asl_pipeline.__version__ import __version__ as v; print('ClinicalASL', v)"

WORKDIR /opt/ClinicalASL/python
ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]
CMD ["/input", "/output", "/workingdir", "--inference-method", "ssvb"]
