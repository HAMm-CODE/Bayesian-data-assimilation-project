# Shared set-up for the Roihu job scripts (sourced, not run on its own).
#
# Python: the Tykky container built for the previous repository, which holds
# PyMC, PyTensor, ArviZ, netCDF4 and a C++ compiler. Our code is NOT installed
# into it (a container cannot be changed); it is found through PYTHONPATH,
# which Python searches before the container's own packages.

export PATH=/projappl/project_2020170/dalec-env/bin:"$PATH"
export PYTHONPATH="$PWD/src${PYTHONPATH:+:$PYTHONPATH}"

# PyTensor keeps compiled code here, so later jobs do not compile again.
export PYTENSOR_FLAGS="base_compiledir=/scratch/project_2020170/pytensor-cache"
mkdir -p /scratch/project_2020170/pytensor-cache jobs/logs results figures

# One thread per chain: the chains already use one CPU each.
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1

echo "=== job ${SLURM_JOB_ID:-local} on $(hostname), $(date) ==="
echo "cpus: ${SLURM_CPUS_PER_TASK:-?}   folder: $PWD"
python - <<'PY'
import sys, importlib
print("python   ", sys.version.split()[0], sys.executable)
for name in ("numpy", "pymc", "pytensor", "arviz", "netCDF4"):
    module = importlib.import_module(name)
    print(f"{name:9s}", getattr(module, "__version__", "?"))
import pytensor, dalec, pathlib
print("C compiler", pytensor.config.cxx or "none (PyTensor would run much slower)")
print("dalec from", pathlib.Path(dalec.__file__).parent)
if not pathlib.Path(dalec.__file__).resolve().is_relative_to(pathlib.Path("src").resolve()):
    sys.exit("dalec is not loaded from this repository's src/ folder; stopping")
PY
echo