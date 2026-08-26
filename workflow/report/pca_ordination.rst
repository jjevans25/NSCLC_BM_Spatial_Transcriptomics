P1-T2 tri-coloured ordination — the figure that motivates every mixed model
downstream (Aim A2). PCA on the top row, UMAP on the bottom, each drawn three
times with the same colourings so the two embeddings can be read against each
other.

**By compartment** and **by site** ask whether the biology separates. **By
patient** asks the question that actually decides model specification: colours
are per subject and each patient's AOIs are joined by a line, so the visual test
is whether a patient's AOIs stay together across compartments rather than which
patient is which. If they do, AOIs within a patient are not independent and the
random intercept ``(1|patient)`` is mandatory rather than stylistic — which is
what P1-T3 then quantifies and P1-T6 records as an ADR.

QC-flagged AOIs are **ringed in red, not removed** (``qc.flag_only``). Whether
they sit apart from their compartment is exactly the question flagging exists to
raise, and it can only be answered by leaving them in.

``TIME-B`` n = 8 throughout.
