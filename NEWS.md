# atmle 0.1.1.9000

* Every fit now includes a variance-floor convex combination of forced- and
  unforced-A HAL candidates for each population target. Its Wald interval uses
  the paired combined influence curve. The default unforced weight cap is 0.5.
* Results identify the estimator and convergence status; both candidates,
  their paired influence curves, and the combined influence curves are saved.
* Targeting avoids dense sample-size-square diagonal matrices, retains explicitly
  forced columns even at zero coefficients, handles zero fluctuation directions,
  and handles scalar singular information matrices.
