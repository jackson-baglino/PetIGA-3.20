"""M3 self-test: planted errors must make the checks fail.

Copies the ledger given by --ledger into a temporary folder, plants one error at a
time (sign, missing factor, wrong exponent, wrong value, dropped premise) in a step,
runs the target check script on the planted copy, and requires the target check to
pass on the original ledger and fail on the planted copy. Runs in parallel.
"""
import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import m3_common as C  # noqa: E402

COVERS = []
HERE = os.path.dirname(os.path.abspath(__file__))

# (name, step, field, old substring, new substring, script, check id that must fail)
PLANTS = [
    ('sign_S9', 'S9', 'statement', '= \\int_a^b', '= -\\int_a^b', 'm3_strong_form.py', 'S9'),
    ('sign_S10', 'S10', 'statement', '= -\\int_a^b', '= \\int_a^b', 'm3_strong_form.py', 'S10'),
    ('sign_S11_source', 'S11', 'statement', '- q\\right)', '+ q\\right)', 'm3_strong_form.py', 'S11'),
    ('sign_S12_source', 'S12', 'statement', '= q \\quad', '= -q \\quad', 'm3_strong_form.py', 'S12'),
    ('factor_S13_storage', 'S13', 'statement', '\\rho c_p \\frac', '\\rho \\frac', 'm3_strong_form.py', 'S13'),
    ('exponent_S14', 'S14', 'statement', '\\frac{\\partial^2 T}{\\partial x^2}', '\\frac{\\partial T}{\\partial x}',
     'm3_strong_form.py', 'S14'),
    ('sign_S15_conduction', 'S15', 'statement', '= k\\,', '= -k\\,', 'm3_strong_form.py', 'S15'),
    ('sign_S15_vs_code', 'S15', 'statement', '= k\\,', '= -k\\,', 'm3_reduction.py', 'reduction_vs_Residual_A1'),
    ('factor_S16', 'S16', 'statement', '\\frac{k}{\\rho c_p}', '\\frac{k}{\\rho}', 'm3_strong_form.py', 'S16'),
    ('factor_S17_source', 'S17', 'statement', '\\frac{q}{\\rho c_p}', '\\frac{q}{\\rho}', 'm3_strong_form.py', 'S17'),
    ('units_S17_source', 'S17', 'statement', '\\frac{q}{\\rho c_p}', '\\frac{q}{\\rho}', 'm3_units.py', 'units_S17'),
    ('grain_in_S15', 'S15', 'statement', '+ q \\quad', '+ q\\,\\psi \\quad', 'm3_reduction.py', 'S18'),
    ('exponent_S20_tau', 'S20', 'statement', '\\frac{\\kappa t}{L^2}', '\\frac{\\kappa t}{L}', 'm3_scaling.py', 'S20'),
    ('factor_S21_kappa', 'S21', 'statement', '\\frac{\\kappa\\,\\Delta T_q}{L^2}', '\\frac{\\Delta T_q}{L^2}',
     'm3_scaling.py', 'S21'),
    ('exponent_S22', 'S22', 'statement', '\\frac{\\Delta T_q}{L^2}\\frac{\\partial^2',
     '\\frac{\\Delta T_q}{L}\\frac{\\partial^2', 'm3_scaling.py', 'S22'),
    ('sign_S23_source', 'S23', 'statement', '+ \\frac{q}{\\rho c_p}', '- \\frac{q}{\\rho c_p}', 'm3_scaling.py', 'S23'),
    ('factor_S24', 'S24', 'statement', '{\\rho c_p\\,\\kappa\\,\\Delta T_q}', '{\\rho c_p\\,\\Delta T_q}',
     'm3_scaling.py', 'S24'),
    ('value_S25_source', 'S25', 'statement', '+ 1 \\quad', '+ 2 \\quad', 'm3_scaling.py', 'S25'),
    ('data_S26', 'S26', 'statement', '\\theta(1,\\tau) = 0', '\\theta(1,\\tau) = 1', 'm3_scaling.py', 'S26'),
    ('exponent_S27_Q', 'S27', 'statement', '\\frac{qL^2}{k\\,T_0}', '\\frac{qL}{k\\,T_0}', 'm3_scaling.py', 'S27'),
    ('factor_S28_source', 'S28', 'statement', '\\frac{q}{\\rho c_p}', '\\frac{q}{\\rho}', 'm3_scaling.py', 'S28'),
    ('factor_S29_Q', 'S29', 'statement', '+ Q,', '+ 2Q,', 'm3_scaling.py', 'S29'),
    ('sign_S30', 'S30', 'statement', '= 1 + Q', '= 1 - Q', 'm3_scaling.py', 'S30'),
    ('value_S31_tc', 'S31', 'statement', 't_c = 1', 't_c = 2', 'm3_scaling.py', 'S31'),
    ('sign_S33_pde', 'S33', 'statement', '= \\kappa\\frac', '= -\\kappa\\frac', 'm3_uniqueness.py', 'S33'),
    ('factor_S34_E', 'S34', 'statement', '\\equiv \\frac12', '\\equiv \\frac14', 'm3_uniqueness.py', 'S34'),
    ('sign_S35', 'S35', 'statement', '= -\\kappa', '= \\kappa', 'm3_uniqueness.py', 'S35'),
    ('sign_S36', 'S36', 'statement', 'T_1 \\equiv T_2', 'T_1 \\equiv -T_2', 'm3_uniqueness.py', 'S36'),
    ('value_S37', 'S37', 'statement', '\\equiv T_0', '\\equiv 2T_0', 'm3_uniqueness.py', 'S37'),
    ('closed_set_S6', 'S6', 'statement', 'T \\in C\\big([0,L]\\times[0,\\infty)\\big)',
     'T \\in C\\big([0,L]\\times(0,\\infty)\\big)', 'm3_uniqueness.py', 'S48'),
    ('closed_set_S6_assump', 'S6', 'statement', 'T \\in C\\big([0,L]\\times[0,\\infty)\\big)',
     'T \\in C\\big([0,L]\\times(0,\\infty)\\big)', 'm3_assumptions.py', 'continuity_at_t0'),
    ('flip_S49', 'S49', 'statement', '(\\xi_*,\\tau_*) \\ge 0', '(\\xi_*,\\tau_*) \\le 0', 'm3_uniqueness.py', 'S49'),
    ('flip_S50', 'S50', 'statement', '\\ge 0', '\\le 0', 'm3_uniqueness.py', 'S50'),
    ('sign_S38', 'S38', 'statement', 'sign of } q', 'sign of } -q', 'm3_uniqueness.py', 'S38'),
    ('drop_S50_from_S38', 'S38', 'from', '"S50"', '"S6"', 'm3_assumptions.py', 'theta_nonnegative'),
    ('drop_S9_from_S11', 'S11', 'from', '"S9", ', '', 'm3_assumptions.py', 'from_lists'),
    ('drop_Tref_const_S1', 'S1', 'statement', 'T_{ref} \\text{ constant in } x \\text{ and } t;', '',
     'm3_assumptions.py', 'tref_constant_for_S13'),
    ('drop_S5_from_S12', 'S12', 'from', ', "S5"', '', 'm3_assumptions.py', 'required_premises'),
    ('drop_S6_from_S21', 'S21', 'from', ', "S6"', '', 'm3_assumptions.py', 'required_premises'),
    ('sign_S39_F0', 'S39', 'statement', 'F_0 = -j(0,t)', 'F_0 = j(0,t)', 'm3_steady_flux.py', 'S39'),
    ('wall_S40', 'S40', 'statement', 'none', 'one term', 'm3_reduction.py', 'S40'),
    ('sign_S41', 'S41', 'statement', '= -1', '= 1', 'm3_steady_flux.py', 'S41'),
    ('factor_S42', 'S42', 'statement', '-\\frac{\\xi^2}{2}', '-\\xi^2', 'm3_steady_flux.py', 'S42'),
    ('value_S43_max', 'S43', 'statement', '\\frac18', '\\frac14', 'm3_steady_flux.py', 'S43'),
    ('factor_S44', 'S44', 'statement', '{2k}', '{k}', 'm3_steady_flux.py', 'S44'),
    ('factor_S44_numeric', 'S44', 'statement', '{2k}', '{k}', 'm3_conservation_energy.py', 'steady_state_vs_S44'),
    ('factor_S45', 'S45', 'statement', '\\frac{qL}{2}', 'qL', 'm3_steady_flux.py', 'S45'),
    ('text_S46', 'S46', 'statement', 'does not apply', 'applies', 'm3_reduction.py', 'S46'),
    ('value_S47', 'S47', 'statement', '\\equiv 0', '\\equiv 1', 'm3_reduction.py', 'S47'),
    ('factor_S51_ic', 'S51', 'statement', 'v(\\xi,0) = \\frac{\\xi(1-\\xi)}{2}', 'v(\\xi,0) = \\xi(1-\\xi)',
     'm3_series.py', 'S51'),
    ('factor_S52', 'S52', 'statement', '\\frac{2((-1)^\\ell - 1)}', '\\frac{((-1)^\\ell - 1)}', 'm3_series.py', 'S52'),
    ('exponent_S53', 'S53', 'statement', '{(\\ell\\pi)^3} = \\begin', '{(\\ell\\pi)^2} = \\begin', 'm3_series.py', 'S53'),
    ('exponent_S54', 'S54', 'statement', '\\frac{4}{(\\ell\\pi)^3}', '\\frac{4}{(\\ell\\pi)^2}', 'm3_series.py', 'S54'),
    ('factor_S55', 'S55', 'statement', '\\frac{4}{\\ell\\pi}', '\\frac{2}{\\ell\\pi}', 'm3_series.py', 'S55'),
    ('value_S56', 'S56', 'statement', '\\ge 0.474', '\\ge 0.48', 'm3_series.py', 'S56'),
    ('exponent_S57_bound', 'S57', 'statement', '\\le \\frac{4}{\\pi}\\exp(-\\pi^2\\tau)',
     '\\le \\frac{4}{\\pi}\\exp(-2\\pi^2\\tau)', 'm3_series.py', 'S57'),
    ('value_S57', 'S57', 'statement', '6.59\\times 10^{-5}', '6.59\\times 10^{-6}', 'm3_series.py', 'S57'),
    ('value_S32', 'S32', 'statement', '6.7\\times 10^{-6}', '6.7\\times 10^{-7}', 'm3_series.py', 'S32'),
    ('value_S32_early', 'S32', 'statement', '\\ge 0.47', '\\ge 0.5', 'm3_series.py', 'S32'),
]


def run(script, ledger):
    p = subprocess.run([sys.executable, os.path.join(HERE, script), '--ledger', ledger],
                       capture_output=True, text=True, timeout=1500)
    res = {}
    for line in p.stdout.splitlines():
        try:
            d = json.loads(line)
        except ValueError:
            continue
        if isinstance(d, dict) and 'check' in d:
            res[d['check']] = d['pass']
    return res


def plant(data, sid, fld, old, new):
    data = json.loads(json.dumps(data))
    step = [s for s in data['steps'] if s['id'] == sid][0]
    if fld == 'from':
        txt = json.dumps(step['from'])
        if old not in txt:
            return None
        step['from'] = json.loads(txt.replace(old, new, 1))
    else:
        if old not in step[fld]:
            return None
        step[fld] = step[fld].replace(old, new, 1)
    return data


def main():
    args = C.cli(__doc__)
    L = C.Ledger(args.ledger)
    rep = C.Reporter('m3_selftest_planted.py', L, {'n_plants': len(PLANTS)})
    tmp = tempfile.mkdtemp(prefix='m3_planted_')
    scripts = sorted({p[5] for p in PLANTS})
    with ThreadPoolExecutor(max_workers=8) as ex:
        base = dict(zip(scripts, ex.map(lambda s: run(s, args.ledger), scripts)))
        jobs = {}
        for name, sid, fld, old, new, script, check in PLANTS:
            data = plant(L.data, sid, fld, old, new)
            if data is None:
                jobs[name] = None
                continue
            path = os.path.join(tmp, name + '.json')
            with open(path, 'w') as f:
                json.dump(data, f)
            jobs[name] = ex.submit(run, script, path)
        planted = {n: (j.result() if j is not None else None) for n, j in jobs.items()}
    for name, sid, fld, old, new, script, check in PLANTS:
        def one(name=name, sid=sid, fld=fld, old=old, new=new, script=script, check=check):
            if planted[name] is None:
                return False, f'pattern {old!r} not found in {sid}.{fld}; cannot plant'
            b = base[script].get(check)
            pl = planted[name].get(check)
            return b is True and pl is False, {'script': script, 'check': check, 'pass_on_original': b,
                                               'pass_on_planted': pl, 'planted': f'{sid}.{fld}: {old!r} -> {new!r}'}
        rep.check('plant_' + name, f'planted error is caught by {script}/{check}', one)
    rep.finish()


if __name__ == '__main__':
    main()
