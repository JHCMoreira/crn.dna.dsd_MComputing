FROM rocker/r-ver:4.4.1

ARG JULIA_VERSION=1.12.5
ARG JULIA_MAJOR_MINOR=1.12

ENV DEBIAN_FRONTEND=noninteractive
ENV JULIA_DEPOT_PATH=/root/.julia
ENV PATH="/opt/julia/bin:${PATH}"

# --------------------------------------------------------------------
# System dependencies
# --------------------------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    gfortran \
    git \
    curl \
    wget \
    ca-certificates \
    pkg-config \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libicu-dev \
    libjpeg-dev \
    libpng-dev \
    libtiff-dev \
    libfribidi-dev \
    libharfbuzz-dev \
    libfreetype6-dev \
    zlib1g-dev \
    python3-pip \
 && rm -rf /var/lib/apt/lists/*

# --------------------------------------------------------------------
# Install Julia
# --------------------------------------------------------------------
RUN curl -fsSL \
    https://julialang-s3.julialang.org/bin/linux/x64/${JULIA_MAJOR_MINOR}/julia-${JULIA_VERSION}-linux-x86_64.tar.gz \
    -o /tmp/julia.tar.gz \
 && mkdir -p /opt \
 && tar -xzf /tmp/julia.tar.gz -C /opt \
 && ln -s /opt/julia-${JULIA_VERSION} /opt/julia \
 && rm /tmp/julia.tar.gz

# --------------------------------------------------------------------
# Install standard R packages (CRAN)
# --------------------------------------------------------------------
RUN install2.r --error --skipinstalled \
    pracma \
    remotes \
    rlist \
    devtools \
    diffeqr \
    deSolve \
    stringr \
    ggplot2 \
    dplyr \
    reshape2 \
    RColorBrewer \
    assertthat \
    magrittr \
    Rdpack \
    gganimate \
    gifski \
    GillespieSSA2 \
    jsonlite \
    testthat

# libSBML is used by the R SBML bridge in R/sbml.R.
RUN python3 -m pip install --no-cache-dir python-libsbml

# --------------------------------------------------------------------
# Install updated JuliaCall from GitHub (Julia 1.12+ compatibility)
# --------------------------------------------------------------------
RUN Rscript -e 'remotes::install_github("JuliaInterop/JuliaCall", upgrade = "never")'

# --------------------------------------------------------------------
# Install & Precompile Julia packages (Added SciMLBase explicitly)
# --------------------------------------------------------------------
RUN julia -e '\
using Pkg; \
Pkg.Registry.update(); \
Pkg.add([ \
    "Suppressor", \
    "RCall", \
    "DiffEqBase", \
    "SciMLBase", \
    "DifferentialEquations" \
]); \
Pkg.precompile();'

# --------------------------------------------------------------------
# Create Julia Startup Hook to auto-load DiffEqBase into Main
# --------------------------------------------------------------------
# This prevents 'UndefVarError: DiffEqBase not defined in Main' in R packages
RUN mkdir -p /root/.julia/config && \
    echo "try; using DiffEqBase, SciMLBase, DifferentialEquations; catch e; @warn \"Could not auto-load SciML packages\" exception=e; end" > /root/.julia/config/startup.jl

# --------------------------------------------------------------------
# Pre-initialize JuliaCall & diffeqr bindings in R
# --------------------------------------------------------------------
RUN Rscript -e 'library(JuliaCall); julia_setup(rebuild = TRUE); library(diffeqr); diffeq_setup()'

CMD ["R"]

# BUILD TIME ~ 30 min
# [+] Building 1898.6s (12/12) FINISHED                                                                                                                         docker:default
#  => [internal] load build definition from Dockerfile                                                                                                                    0.0s
#  => => transferring dockerfile: 3.19kB                                                                                                                                  0.0s
#  => [internal] load metadata for docker.io/rocker/r-ver:4.4.1                                                                                                           1.3s
#  => [internal] load .dockerignore                                                                                                                                       0.2s
#  => => transferring context: 207B                                                                                                                                       0.1s
#  => CACHED [1/8] FROM docker.io/rocker/r-ver:4.4.1@sha256:78cb94ce2db23aaaf7b546450fcf70b5a3f2ace5a9b5fa1f87217da329211312                                              0.1s
#  => => resolve docker.io/rocker/r-ver:4.4.1@sha256:78cb94ce2db23aaaf7b546450fcf70b5a3f2ace5a9b5fa1f87217da329211312                                                     0.1s
#  => [2/8] RUN apt-get update && apt-get install -y --no-install-recommends     build-essential     gfortran     git     curl     wget     ca-certificates     pkg-con  67.1s
#  => [3/8] RUN curl -fsSL     https://julialang-s3.julialang.org/bin/linux/x64/1.12/julia-1.12.5-linux-x86_64.tar.gz     -o /tmp/julia.tar.gz  && mkdir -p /opt  && ta  45.4s 
#  => [4/8] RUN install2.r --error --skipinstalled     pracma     remotes     rlist     devtools     diffeqr     deSolve     stringr     ggplot2     reshape2     RCol  308.1s 
#  => [5/8] RUN Rscript -e 'remotes::install_github("JuliaInterop/JuliaCall", upgrade = "never")'                                                                        25.9s 
#  => [6/8] RUN julia -e 'using Pkg; Pkg.Registry.update(); Pkg.add([     "Suppressor",     "RCall",     "DiffEqBase",     "SciMLBase",     "DifferentialEquations" ]  1252.3s 
#  => [7/8] RUN mkdir -p /root/.julia/config &&     echo "try; using DiffEqBase, SciMLBase, DifferentialEquations; catch e; @warn "Could not auto-load SciML packages" e  0.5s 
#  => [8/8] RUN Rscript -e 'library(JuliaCall); julia_setup(rebuild = TRUE); library(diffeqr); diffeq_setup()'                                                           27.7s 
#  => exporting to image                                                                                                                                                169.4s 
#  => => exporting layers                                                                                                                                               104.8s 
#  => => exporting manifest sha256:47a8e4eff1c024d434d60f796eb7dfd4e2e04090de5b98e49dd86fabb8609c93                                                                       0.1s 
#  => => exporting config sha256:354cc22e40ae6c2cca9f1fcb64f2dbaa1f2030d649d149128314d25f743748b4                                                                         0.1s 
#  => => exporting attestation manifest sha256:412f8ad5c4b8039e9d48f2bb4e9fa338c43363636230c7591d3a6ceefd48b277                                                           0.2s 
#  => => exporting manifest list sha256:28a6b0b529a2a0731a23500cf008316e43a30f72976c678db24acaedcc92b12b                                                                  0.1s
#  => => naming to docker.io/library/dnar-dev:latest                                                                                                                      0.0s
#  => => unpacking to docker.io/library/dnar-dev:latest                                                                                                                  63.9s
