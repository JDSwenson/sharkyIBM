# sharkyIBM Development Documentation

This `docs/` folder contains development notes, memory files, and conversation exports related to the sharkyIBM package. These files are **excluded from package builds** and are intended for developers and maintainers only.

## Contents

- **AGENTS.md** — Project memory: design decisions, function specifications, test results, and implementation notes. Updated after major development milestones.
- **CITATION_GUIDE.md** — Guide for maintaining package citations and BibTeX bibliography.
- **DOCUMENTATION_ROADMAP.md** — Map of parameter documentation across the three main functions.
- **DOCUMENTATION_COMPLETION_LOG.md** — Detailed account of roxygen2 documentation work and quality assurance.
- **SPECIES_GENERALIZATION_TODO.md** — Roadmap for expanding the package to other species beyond dolphins.
- **conversation-export.{html,md}** — Exports from development conversations (large, non-essential files).

## Quick Start for Maintainers

If you're picking up this project:

1. Read **AGENTS.md** first to understand the design, feature status, and last known state.
2. Check **DOCUMENTATION_ROADMAP.md** to locate parameter documentation.
3. Refer to **CITATION_GUIDE.md** when updating references.
4. See **SPECIES_GENERALIZATION_TODO.md** for planned extensions.

## Updating Memory Files

After significant changes or feature additions:

```bash
# You can update AGENTS.md via the `/savememory` command in Posit Assistant
```

Or manually edit the file to record:
- What was changed
- Why it was changed
- Test results
- Next steps

---

**Last Updated:** September 2, 2026
