# Two targets:
#   docker build --target verify -t cap-verify . && docker run --rm cap-verify   # level 1, standard library only
#   docker build --target full   -t cap-full   . && docker run --rm cap-full     # level 2: build evans, recompute samples
# Level 3 (everything, ~35 CPU-hours): docker run --rm cap-full sh scripts/full_recompute.sh g0 consts pi g1 thmA lamstar lambeta g2

# ---------------------------------------------------------------- level 1
FROM python:3.12-slim AS verify
WORKDIR /cap
COPY . /cap
CMD ["sh", "scripts/quick_verify.sh"]

# ---------------------------------------------------------------- levels 2 and 3
# Ubuntu 24.04 ships FLINT 3.0.1, the version used for evans; the package version is pinned.
# If that exact package is ever unavailable, build FLINT 3.0.1 from source instead:
#   https://github.com/flintlib/flint/releases/download/v3.0.1/flint-3.0.1.tar.gz
#   ./configure --prefix=/usr/local && make -j && make install   (then: make FLINT_PREFIX=/usr/local in g2_evans)
FROM ubuntu:24.04 AS full
ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential make libflint-dev=3.0.1-3.1build1 libgmp-dev libmpfr-dev \
        python3 python3-pip python3-venv ca-certificates \
    && rm -rf /var/lib/apt/lists/*
RUN python3 -m venv /opt/venv && /opt/venv/bin/pip install --no-cache-dir \
        python-flint==0.9.0 numpy==2.5.3 scipy==1.18.1 mpmath==1.3.0
ENV PATH=/opt/venv/bin:$PATH
WORKDIR /cap
COPY . /cap
# same flags as the original build (-O2 -march=native); the sample comparison of subdivision counts assumes the same
# floating-point behaviour (x86-64 with FMA). Rigour does not depend on it.
RUN make -C g2_evans evans
CMD ["sh", "-c", "sh scripts/quick_verify.sh && python3 scripts/recompute_samples.py"]
