from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Tuple, Optional, Dict
import argparse


# -------------------------
# 1) Configuration “ctx”
# -------------------------
@dataclass
class SimulationConfig:
    # Simulation parameters
    dim: int = 2
    temp0: float = 273.15 - 30  # K
    grad_temp0: Tuple[float, float, float] = (10.0, 0.0, 0.0)  # K/m in x, y, z

    init_mode: str = "layered"  # "layered", "circle", or anything else to read from file

    # Required at runtime (populated from CLI and .env)
    init_dir: Path | None = None
    Nx: Optional[int] = None
    Ny: Optional[int] = None
    Nz: Optional[int] = None
    Lx: Optional[float] = None
    Ly: Optional[float] = None
    Lz: Optional[float] = None
    eps: Optional[float] = None


# -------------------------
# 2) .env discovery & parsing
# -------------------------
ENV_CANDIDATE_NAMES = ("grains.env",)


def find_env_file(init_dir: Path) -> Path:
    """Locate a grains.env file inside init_dir.

    Priority:
      1) A file literally named "grains.env"
      2) If not present, exactly one file matching "*.env"
    """
    if not init_dir.is_dir():
        raise SystemExit(f"❌ init_dir is not a directory: {init_dir}")

    # 1) Prefer a literal "grains.env"
    for name in ENV_CANDIDATE_NAMES:
        candidate = init_dir / name
        if candidate.is_file():
            return candidate

    # 2) Otherwise, look for a single *.env
    envs = sorted(init_dir.glob("*.env"))
    if len(envs) == 1:
        return envs[0]
    elif len(envs) == 0:
        raise SystemExit(f"❌ No grains.env found in {init_dir}. Expected a file named 'grains.env' or exactly one '*.env'.")
    else:
        listing = "\n  - " + "\n  - ".join(str(p) for p in envs)
        raise SystemExit(f"❌ Multiple .env candidates found in {init_dir}:{listing}\nSettle on a single file named 'grains.env'.")


def parse_simple_env(path: Path) -> Dict[str, str]:
    """Minimal .env parser: KEY=VALUE per line, ignores blanks and # comments."""
    kv: Dict[str, str] = {}
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith('#'):
            continue
        # Allow inline comments (KEY=VALUE  # comment)
        if '#' in line:
            line = line.split('#', 1)[0].strip()
        if '=' not in line:
            continue
        k, v = line.split('=', 1)
        k = k.strip()
        v = v.strip().strip('"\'')
        if k:
            kv[k] = v
    return kv


def coerce_mesh_params(kv: Dict[str, str]) -> Dict[str, object]:
    """Coerce required params to proper types, raise with clear error if missing."""
    required = ["Nx", "Ny", "Nz", "Lx", "Ly", "Lz", "eps"]
    missing = [k for k in required if k not in kv]
    if missing:
        raise SystemExit(f"❌ Missing keys in .env: {', '.join(missing)}")

    try:
        return {
            "Nx": int(kv["Nx"]),
            "Ny": int(kv["Ny"]),
            "Nz": int(kv["Nz"]),
            "Lx": float(kv["Lx"]),
            "Ly": float(kv["Ly"]),
            "Lz": float(kv["Lz"]),
            "eps": float(kv["eps"]),
        }
    except ValueError as e:
        raise SystemExit(f"❌ Type conversion error while reading .env: {e}")


# -------------------------
# 3) CLI wiring
# -------------------------

def parse_args() -> SimulationConfig:
    parser = argparse.ArgumentParser(description="Simulation configuration")
    parser.add_argument("init_dir", type=Path, help="Directory that contains the .env with Nx,Ny,Nz,Lx,Ly,Lz,eps")
    args = parser.parse_args()

    cfg = SimulationConfig(init_dir=args.init_dir)

    env_path = find_env_file(cfg.init_dir)
    kv = parse_simple_env(env_path)
    coerced = coerce_mesh_params(kv)

    cfg.Nx = coerced["Nx"]  # type: ignore[assignment]
    cfg.Ny = coerced["Ny"]  # type: ignore[assignment]
    cfg.Nz = coerced["Nz"]  # type: ignore[assignment]
    cfg.Lx = coerced["Lx"]  # type: ignore[assignment]
    cfg.Ly = coerced["Ly"]  # type: ignore[assignment]
    cfg.Lz = coerced["Lz"]  # type: ignore[assignment]
    cfg.eps = coerced["eps"]  # type: ignore[assignment]

    return cfg


if __name__ == "__main__":
    config = parse_args()
    # Quick echo so you can verify what was read
    print("init_dir:", config.init_dir)
    print("Nx,Ny,Nz:", config.Nx, config.Ny, config.Nz)
    print("Lx,Ly,Lz,eps:", config.Lx, config.Ly, config.Lz, config.eps)