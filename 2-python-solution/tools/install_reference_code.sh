#!/usr/bin/env bash
# Installs the existing NICER light-curve code X-PSI (the reference for the stopwatch) into ONE folder on the Linux side:
# /root/physics_nicer. It puts its own Python (Miniforge) inside that folder. Nothing outside the folder is changed:
# no system package is installed (if the system has no C compiler, one is put inside X-PSI's own environment).
# To undo: delete /root/physics_nicer.
# Follows X-PSI's documented steps (docs/source/install.rst): conda env create -f environment.yml ; CC=gcc pip install .
# NOT RUN AND NOT TESTED when written (22 Nov 2026).
set -u
BASE=/root/physics_nicer
mkdir -p "$BASE" && cd "$BASE" || { echo "INSTALL_FAILED cannot make $BASE"; exit 10; }
log() { echo "[$(date -u +%H:%M:%S)] $*"; }
{
log "start"
if [ ! -x "$BASE/miniforge/bin/conda" ]; then
  log "getting Miniforge (its own Python, inside this folder)"
  curl -sSL -m 900 -o miniforge.sh https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh || { echo "INSTALL_FAILED download"; exit 11; }
  bash miniforge.sh -b -p "$BASE/miniforge" || { echo "INSTALL_FAILED miniforge"; exit 12; }
fi
export PATH="$BASE/miniforge/bin:$PATH"
if [ ! -d "$BASE/xpsi" ]; then
  log "getting X-PSI's source"
  git clone https://github.com/xpsi-group/xpsi.git "$BASE/xpsi" || { echo "INSTALL_FAILED clone"; exit 13; }
  # the same revision whose source was read for the specification
  ( cd "$BASE/xpsi" && git checkout -q bb5f461e605b22c6251b155b27a914a558757f93 ) || log "note: could not check out the revision that was read; using the newest"
fi
if [ ! -d "$BASE/miniforge/envs/xpsi" ]; then
  log "making X-PSI's environment the documented way (the long step)"
  if ! conda env create -q -f "$BASE/xpsi/environment.yml"; then
    log "the documented way failed; making the environment from the conda-forge channel only, same packages"
    conda env remove -q -y -n xpsi >/dev/null 2>&1
    conda create -q -y -n xpsi -c conda-forge "python>=3.10" numpy "cython~=3.0.11" matplotlib scipy pandas wrapt "setuptools>=61" gsl pytest getdist tqdm h5py nestcheck astropy emcee pip \
      || { echo "INSTALL_FAILED environment"; exit 14; }
  fi
fi
if command -v gcc >/dev/null 2>&1; then
  XCC=gcc
  log "using the system's C compiler: $(gcc --version | head -n 1)"
else
  log "no C compiler on the system; putting one inside X-PSI's environment"
  conda install -q -y -n xpsi -c conda-forge gcc_linux-64 || { echo "INSTALL_FAILED compiler"; exit 17; }
  XCC=x86_64-conda-linux-gnu-gcc
fi
log "building X-PSI with CC=$XCC"
( cd "$BASE/xpsi" && conda run -n xpsi env CC="$XCC" pip install -q . ) || { echo "INSTALL_FAILED build"; exit 15; }
mkdir -p "$BASE/run" && cd "$BASE/run"
conda run -n xpsi python -c "import xpsi; print('xpsi', xpsi.__version__)" || { echo "INSTALL_FAILED import"; exit 16; }
du -sh "$BASE" | sed 's/^/size: /'
log "end"
echo INSTALL_DONE
} 2>&1 | tee -a "$BASE/install.log"
