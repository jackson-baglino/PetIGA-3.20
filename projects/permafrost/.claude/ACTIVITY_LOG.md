
---

**Session ended:** 2026-04-17 14:02:23


---

**Session ended:** 2026-04-17 11:55:47


## 2026-04-16 — Post-processing updates and bug fixes

- Updated `postprocess/analyze_interface.py`: read sediment DOF (sol[:,3]), fix air = 1−φ_i−φ_s, add sediment volume fraction to metrics and 4×2 panel plot
- Updated `postprocess/plot2D_snapshot.py`: expand `_plot_cuts` from 2×3 to 2×4 to include sediment cross-sections
- Updated `postprocess/plotpermafrost.py`: add SedPhase DOF export, fix file-exists path bug, add CLI flags
- Updated `scripts/plotpermafrost.py`: fix `❌ Error processing` bug caused by erroneous `np.newaxis` in air-phase computation
- Set up `.claude/ACTIVITY_LOG.md` and Stop hook for session logging; added activity log instruction to `CLAUDE.md`

---

**Session ended:** 2026-04-16 21:55:05


---

**Session ended:** 2026-04-16 21:51:32


---

**Session ended:** 2026-04-16 21:45:52


---

**Session ended:** 2026-04-16 21:42:34

# Activity Log — Permafrost Project

Newest entries appear at the top. Each session is separated by `---`.

---

## 2026-04-16 — Configure activity log

Set up `.claude/ACTIVITY_LOG.md` and a Stop hook to automatically timestamp each session end. Added instruction to `CLAUDE.md` requiring a summary entry before stopping.

---