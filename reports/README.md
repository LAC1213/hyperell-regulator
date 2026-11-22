# Beilinson survey: conductor < 1000, 100 digits

- [Saved data, exact divisor coefficients and interval endpoints](beilinson_100_digits.json).
- [Tables of divisors, packet torsion orders, endomorphism types and numerical comparisons](beilinson_100_digits.tex).

The snapshot contains 87 records and 41 completed regulators.
Failed attempts remain in the JSON. Q and CM refer to rational endomorphism algebra labels.
Divisors use the original LMFDB coordinates; the integration charts are recorded separately. Packet order is the exponent of point differences,
not packet cardinality. Unknown exact orders are left blank rather than inferred from Jacobian torsion.

Regulator determinants use weights `(1, 2 Re(z), |z|^2)` and full cycle periods.
Certified interval endpoints and their analytic assumptions are retained; numerical agreement does not prove equality.

Regenerate from a completed experiment without repeating integrations:

```sh
sage -python -m hyperell_regulator.experiments.beilinson_report --input results/beilinson_results_certified.json --bound 1000 --digits 100
```

Source SHA-256: `37e3fedbf704d45cdd920707f77ad3c1746099ffb6e47a6c25e6d91a477fa682`.
