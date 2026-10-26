# When does the early, width-dependent transient end? (2026-10-25)

`relaxation_end.py` reads the paired interface-width runs (phi 0.325, -20 C,
seeds 1701-1703, eps = 1 and 2 um) from the campaign directory and computes,
for each baseline time t_b, the rise k_iso(30 d)/k_iso(t_b) - 1 at both widths.
Output: `relaxation_end.txt`.

The two widths disagree by 9.4 points at t_b = 1 d (theta 11 at 1 um), 1.8 at
2.7 d (theta 30), 1.0 at 4 d, 0.6 at 5 d and 0.0 at 8 d. Seed scatter of the
rise is about 2 points. "11 tau_sub" was never a measured end of relaxation:
it is one day at -20 C written in tau_sub (introduced 2026-09-27 as a common
baseline across temperatures). Only eps = 2 um is tested against eps = 1 um,
so where the 1 um run itself converges is not known: theta about 22 if the
transient scales with tau_sub, later if it does not.
