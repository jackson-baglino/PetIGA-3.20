# ekcfg.py — tiny, fast config for ek_eff.py (FEniCSx)
from __future__ import annotations
import os, json
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Tuple, Optional, Mapping

# ---------- small parsers (no external deps) ----------
_TRUE = {"1", "true", "t", "yes", "y", "on"}
_FALSE = {"0", "false", "f", "no", "n", "off"}

def _as_bool(s: str, default: bool=False) -> bool:
    if s is None: return default
    s2 = s.strip().lower()
    if s2 in _TRUE: return True
    if s2 in _FALSE: return False
    raise ValueError(f"Bad boolean env value: {s!r}")

def _as_int(s: str, default: int) -> int:
    return int(s) if s not in (None, "") else default

def _as_float(s: str, default: float) -> float:
    return float(s) if s not in (None, "") else default

def _as_str(s: str, default: str) -> str:
    return s if s not in (None, "") else default

def _as_vec3(s: str, default: Tuple[float,float,float]) -> Tuple[float,float,float]:
    if s in (None, ""): return default
    parts = [p.strip() for p in s.replace(";", ",").split(",") if p.strip()!=""]
    if len(parts) != 3: raise ValueError(f"Expected 3 components, got {parts}")
    return (float(parts[0]), float(parts[1]), float(parts[2]))

# ---------- the config ----------
@dataclass(frozen=True)
class Config:
    # Domain & mesh
    dim: int              = 2
    Lx: float             = 1.0e-3
    Ly: float             = 1.0e-3
    Lz: float             = 1.0e-3
    Nx: int               = 512
    Ny: int               = 512
    Nz: int               = 1
    p:  int               = 1              # polynomial degree (CG)

    # Material
    thcond_ice: float     = 2.29
    thcond_air: float     = 0.02
    eps: float            = 0.01

    # Macroscopic fields (only if you keep them)
    temp0: float          = 0.0
    grad_temp0: Tuple[float,float,float] = (0.0, 0.0, 0.0)

    # Forcing/BC parameters you may still want to keep for compatibility
    T_top: float          = 0.0
    q_bottom: float       = 0.0

    # Initialization
    init_mode: str        = "file"         # "circle" | "layered" | "file"
    init_dir: Path        = Path(".")
    sol_index: int        = -1

    # Output
    outdir: Path          = Path("./k_eff_output")
    output_binary: bool   = False

    # ---------- helpers / derived ----------
    @property
    def vol(self) -> float:
        return (self.Lx * self.Ly) if self.dim == 2 else (self.Lx * self.Ly * self.Lz)

    @property
    def mesh_shape(self) -> Tuple[int,int,int]:
        return (self.Nx, self.Ny, (1 if self.dim==2 else self.Nz))

    def with_overrides(self, **kw) -> "Config":
        """Create a copy with specified fields changed."""
        return replace(self, **kw)

    def to_json(self, path: Path) -> None:
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w") as f:
            json.dump({
                **{k: getattr(self, k) for k in self.__dataclass_fields__.keys()},
                "init_dir": str(self.init_dir),
                "outdir":   str(self.outdir),
            }, f, indent=2)

    @staticmethod
    def from_json(path: Path) -> "Config":
        data = json.loads(path.read_text())
        data["init_dir"] = Path(data["init_dir"])
        data["outdir"]   = Path(data["outdir"])
        return Config(**data)

    def summary(self) -> str:
        ms = f"{self.Nx}x{self.Ny}" + ("" if self.dim==2 else f"x{self.Nz}")
        return (f"[dim={self.dim}] L=({self.Lx},{self.Ly}" +
                ("" if self.dim==2 else f",{self.Lz}") +
                f") N={ms} p={self.p} eps={self.eps} mode={self.init_mode} idx={self.sol_index}")

    # ---------- construction from environment ----------
    @classmethod
    def from_env(cls, env: Optional[Mapping[str,str]]=None, base: Optional["Config"]=None) -> "Config":
        """Build from os.environ (or provided mapping). Bash-exported vars override defaults."""
        env = os.environ if env is None else env
        cfg = base or cls()

        try:
            cfg = cfg.with_overrides(
                dim = _as_int(env.get("dim"), cfg.dim),
                Lx  = _as_float(env.get("Lx"), cfg.Lx),
                Ly  = _as_float(env.get("Ly"), cfg.Ly),
                Lz  = _as_float(env.get("Lz"), cfg.Lz),

                Nx  = _as_int(env.get("Nx"), cfg.Nx),
                Ny  = _as_int(env.get("Ny"), cfg.Ny),
                Nz  = _as_int(env.get("Nz"), cfg.Nz),

                p   = _as_int(env.get("P"), cfg.p),  # optional

                thcond_ice = _as_float(env.get("THCOND_ICE"), cfg.thcond_ice),
                thcond_air = _as_float(env.get("THCOND_AIR"), cfg.thcond_air),
                eps        = _as_float(env.get("eps"), cfg.eps),

                temp0      = _as_float(env.get("TEMP0"), cfg.temp0),
                grad_temp0 = _as_vec3(env.get("GRAD_TEMP0"), cfg.grad_temp0),

                T_top      = _as_float(env.get("TEMP_TOP"), cfg.T_top),
                q_bottom   = _as_float(env.get("FLUX_BOTTOM"), cfg.q_bottom),

                init_mode  = _as_str(env.get("INIT_MODE"), cfg.init_mode),
                init_dir   = Path(_as_str(env.get("INIT_DIR"), str(cfg.init_dir))),

                sol_index  = _as_int(env.get("SOL_INDEX"), cfg.sol_index),

                outdir         = Path(_as_str(env.get("OUTPUT_DIR"), str(cfg.outdir))),
                output_binary  = _as_bool(env.get("OUTPUT_BINARY"), cfg.output_binary),
            )
        except Exception as e:
            raise RuntimeError(f"Invalid environment configuration: {e}") from e

        return cfg._validated()

    # ---------- validation & normalization ----------
    def _validated(self) -> "Config":
        if self.dim not in (2, 3):
            raise ValueError("dim must be 2 or 3")
        if self.dim == 2:
            # normalize 2D-friendly values
            nz = 1
            lz = 1.0 if self.Lz == 0 else self.Lz
            return replace(self, Nz=nz, Lz=lz)
        return self
    
def thermal_cond(cfg, ice: float) -> tuple[float, float]:
    """
    Compute effective thermal conductivity and its derivative wrt ice content.

    Parameters
    ----------
    cfg : Config or AppCtx
        Holds material properties `thcond_ice` and `thcond_air`.
    ice : float
        Ice volume fraction (can be outside [0,1]; clamped internally).

    Returns
    -------
    cond : float
        Effective conductivity [W/m·K].
    dcond_ice : float
        Derivative of conductivity with respect to ice fraction.
    """
    dice = 1.0
    dair = 1.0
    air = 1.0 - ice

    # Clamp ice and air to [0,1], update derivative flags
    if ice < 0.0:
        ice = 0.0
        dice = 0.0
    if air < 0.0:
        air = 0.0
        dair = 0.0

    cond_ice = cfg.thcond_ice
    cond_air = cfg.thcond_air

    cond = ice * cond_ice + air * cond_air
    dcond_ice = cond_ice * dice - cond_air * dair

    return cond, dcond_ice

