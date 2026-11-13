# Further simplifications — implemented after approval

Both proposed simplifications were explicitly approved and have been applied.

1. **One full-period quadrature engine.** Float64 cycle periods and signed
   minors now delegate to the guarded engine in `regulator/precision.py`,
   converting only the final results to NumPy arrays or Python floats.
   The separate SciPy cycle integration and duplicate minor assembly have
   been removed. Each MST cycle integrates one branch-to-branch edge, so
   the old multi-segment accumulation and corner-continuation state are
   also unnecessary. The independent two-dimensional plane integrator
   remains available for reference calculations.

   The existing `plane.period_matrix`, `cplus` and `cminus` entry points
   retain their tolerance arguments. Positive tolerances select an internal
   decimal target of at least 15 digits; nonpositive values disable that
   tolerance. At least one finite positive tolerance is required. Return
   types remain complex128 NumPy arrays or Python floats.

2. **Remove unused experiment API arguments.** Numerical column helpers
   no longer accept `order`, and `beilinson_determinant_stokes` now takes
   `(G, imP, imQ, *, digits=None, pole="infinity")` without `nP` or `nQ`.
   Triple-ranking helpers no longer accept `mp`. Callers in the conductor
   and family experiments, tests and documentation use the new signatures.
   Numerical options are keyword-only to prevent old order arguments from
   silently binding to precision or pole options.

   Torsion orders remain in exact certification, component checks and report
   metadata. The MST cycle basis and mathematical definitions are retained.
