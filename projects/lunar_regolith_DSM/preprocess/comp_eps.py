"""comp_eps.py: thin wrapper kept so the old command line keeps working.

The logic now lives in the agent-facing tool  agentic/tools/mesh_params/  (typed input and output, structured checks,
tests that pin it to the previous numbers). Same flags as before, plus --json / --json-file / --spec / --as-json:

  python preprocess/comp_eps.py --Lx 6e-4 --Ly 6e-4 --Lz 6e-4 --T0 -10 --alpha 1e-3
  python preprocess/comp_eps.py --Lx 6e-4 --T0 -20 --alpha_range 1e-4 1e-1 12
  python preprocess/comp_eps.py --Lx 6e-4 --T0 -5 --alpha 1e-3 --patch inputs/geometry/<family>/<name>.opts   # sets -eps, -Nx/-Ny/-Nz only
  python preprocess/comp_eps.py --Lx 6e-4 --T0 -5 --alpha 1e-3 --quiet

--patch no longer writes -M_sub, -alpha_sub, -tau_sub, -lam_sub, -beta_sub or -d0_sub: the solver does not read those options.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))   # project root, so `agentic` imports

from agentic.tools.mesh_params.cli import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())
