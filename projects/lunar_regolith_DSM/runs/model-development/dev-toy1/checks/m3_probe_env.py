"""M3 check: the Python environment has the packages the M3 checks use.

Prints one JSON line; exits 0 only if sympy, numpy, scipy and pint import.
(No LaTeX parser backend is needed: m3_common.py has its own parser.)
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m3_common as C  # noqa: E402

COVERS = []


def main():
    args = C.cli(__doc__)
    L = C.Ledger(args.ledger)
    rep = C.Reporter('m3_probe_env.py', L, {'python': sys.version.split()[0]})

    def env():
        out = {}
        for mod in ('sympy', 'numpy', 'scipy', 'pint'):
            try:
                m = __import__(mod)
                out[mod] = getattr(m, '__version__', 'present')
            except Exception as exc:  # noqa: BLE001
                out[mod] = f'missing: {exc}'
        return all(not str(v).startswith('missing') for v in out.values()), out
    rep.check('environment', 'sympy, numpy, scipy and pint are importable', env)
    rep.finish()


if __name__ == '__main__':
    json  # noqa: B018 - keep import for clarity of output format
    main()
