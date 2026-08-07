# inputs/scratch — quarantine, not an archive

Every geometry file, experiment file, packing and mesh that is not currently in
use, kept in its original sub-directory layout so it can be moved back
unchanged.

**Before writing a new input file, search here first:**

```bash
find inputs/scratch -name '*<something>*'
```

If it exists, pull it back out rather than writing an equivalent under a new
name:

```bash
git mv inputs/scratch/geometry/<family>/<name>.opts inputs/geometry/<family>/<name>.opts
```

The run scripts do this search for you — when a geometry or experiment name
misses, they print the exact `git mv` to run.

**This directory will be cleared out wholesale.** Anything worth keeping has to
be back in the live tree before then.
