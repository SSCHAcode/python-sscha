# SSCHA Docker Images

Two images are provided:

- **`sscha-cuda`** — base image with SSCHA, PyTorch, MPI, and Julia
- **`sscha-freeenergy`** — extends `sscha-cuda` with `nequip` and `openequivariance` for the free energy minimization workflow

## Requirements

- Docker with [nvidia-container-toolkit](https://docs.nvidia.com/datacenter/cloud-native/container-toolkit/install-guide.html)
- NVIDIA driver **>= 560** (CUDA 12.6 support)

Verify your driver version:
```bash
nvidia-smi
```

## Building

The images must be built in order since `sscha-freeenergy` depends on `sscha-cuda`. Both commands are run from the **repository root**.

```bash
docker build -f docker/Dockerfile.conda -t sscha-cuda .
docker build -f docker/Dockerfile.freeenergy -t sscha-freeenergy .
```

Use `--no-cache` to force a clean rebuild:
```bash
docker build --no-cache -f docker/Dockerfile.conda -t sscha-cuda .
docker build -f docker/Dockerfile.freeenergy -t sscha-freeenergy .
```

## Setup

Before running the example, obtain the `pdh_kitware` folder and place it in the repository root. The file structure should then resemble the following:

```
python-sscha/  <-- Repo root
└── pdh_kitware/
    └── sscha/
        ├── run_sscha.py
        └── 4x4x4_harmonic_dyn*
```

## Running the PdH free energy minimization example

The example is in `pdh_kitware/sscha/`. It requires the NequIP-OAM-XL model, which is downloaded and compiled automatically on first launch.

Mount a named volume for the model so it persists across runs and is only compiled once:

```bash
docker run --rm --gpus all \
    -v nequip-model:/opt/nequip \
    -v "$(pwd)/pdh_kitware/sscha":/workspace \
    -w /workspace \
    sscha-freeenergy -c "python run_sscha.py"
```

On first launch you will see:
```
Compiling NequIP-OAM-XL model (first run only)...
Model compiled to /opt/nequip/NequIP-OAM-XL.nequip.pt2
```

Subsequent runs skip compilation and go straight to the calculation. Output files (`data/`, `minim_info*`) are written into `pdh_kitware/sscha/` on the host.

### Using a different model path

Override `NEQUIP_MODEL_PATH` to use a pre-compiled model stored elsewhere:

```bash
docker run --rm --gpus all \
    -v /path/to/models:/models \
    -v "$(pwd)/pdh_kitware/sscha":/workspace \
    -w /workspace \
    -e NEQUIP_MODEL_PATH=/models/NequIP-OAM-XL.nequip.pt2 \
    sscha-freeenergy -c "python run_sscha.py"
```

## Interactive shell

```bash
docker run -it --gpus all \
    -v nequip-model:/opt/nequip \
    -v "$(pwd)/pdh_kitware/sscha":/workspace \
    -w /workspace \
    sscha-freeenergy
```

## Verifying the build

```bash
# SSCHA imports
docker run --rm sscha-cuda -c "python -c 'import sscha; import cellconstructor; print(\"SSCHA OK\")'"

# CUDA visible to PyTorch
docker run --rm --gpus all sscha-cuda -c "python -c 'import torch; print(torch.cuda.is_available())'"

# MPI
docker run --rm sscha-cuda -c "mpirun -n 2 python -c 'from mpi4py import MPI; print(MPI.COMM_WORLD.Get_rank())'"

# nequip + openequivariance
docker run --rm --gpus all sscha-freeenergy -c "python -c 'import nequip; import openequivariance; print(\"OK\")'"
```
