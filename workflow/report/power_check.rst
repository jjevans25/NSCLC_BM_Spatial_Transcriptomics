P0-T8 power reality-check for the project's central contrast: the tumour immune
microenvironment in the primary lung tumour (``TIME-L``, n = 15) against the
brain metastasis (``TIME-B``, **n = 8**).

Curves are simulation-based, using the design exactly as ``samples.tsv`` records
it rather than an idealised balanced version — 16 distinct patients, two of whom
contribute two ``TIME-L`` AOIs, and **five of whom appear in both groups**. Those
five induce a within-subject correlation across the groups being compared, which
a two-sample t-test on 15 vs 8 would assume away and so report a more optimistic
number than the design supports.

Effects are standardised so that patient and residual variance sum to 1; the
x-axis therefore reads in SD units and ICC is the share of variance between
patients. Significance is judged against a t distribution with 14 degrees of
freedom, a conservative stand-in for the Satterthwaite approximation ``lmerTest``
will use in Phase 2; the Wald z that ``MixedLM`` reports natively is recorded in
``results/tables/power.tsv`` alongside, and is meaningfully more optimistic at
this sample size.

The minimum detectable effect at 80% power is annotated per ICC in the legend.
Read it as the floor on what Phase 2 can claim about brain immune contexture.
