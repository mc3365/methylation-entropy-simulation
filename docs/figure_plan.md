# Figure plan for the per-CpG conditional null

Do not draw a single stochastic entropy band against pooled mean methylation.
Windows with the same pooled mean can have different four-CpG methylation
profiles, read depths, and therefore different null distributions.

## Small tool figure

Show each four-CpG window's observed entropy relative to its own null using
`entropy_minus_independent_null` or its null percentile. Plot the effect against
mean methylation and show read depth by facets or a separate depth panel.
Distinguish a low observed entropy (coordination) from a high observed entropy
(dispersion). A few illustrative windows may include their individual null
histograms, exact CpG rates, observed value, and empirical tail probability.
Explain that each null permutes the four CpGs independently across the same
number of reads while fixing the count at every CpG.

## DMV biological figure

Aggregate window effects within each DMV and biological replicate, with a
minimum informative-window rule and coverage reported alongside the effect.
Plot replicate-level DMV summaries across unpassaged, P0, and P5. Use a
predefined DMV set and comparison rule. Do not treat overlapping windows as
independent replicates, and do not label a low-entropy result as proof of a
specific regulatory mechanism.

Generate these figures only after the real-data pilot establishes the usable
coverage range and the chosen DMV aggregation rule. Record excluded windows
and show uncertainty at the biological-replicate level.
