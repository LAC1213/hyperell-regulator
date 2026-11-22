Curve explorer
==============

`Open the explorer <explorer/index.html>`_ to inspect the Stokes computation
for the selected degree-zero divisors on a genus two curve
:math:`C:y^2+h(x)y=f(x)` in its original LMFDB coordinates. The catalog contains the 41 successful comparisons
from the recent 100-digit survey below conductor 1000, including five
curves labelled ``CM`` in the rational endomorphism algebra field.

The viewer
----------

Choose a curve to see its equation, selected divisors, packet torsion order
and rational endomorphism algebra label. Packet order is the exponent
:math:`N=\operatorname{lcm}_P\operatorname{ord}(P-P_0)` in the Jacobian, rather than
the number of points or the order of the ambient Jacobian torsion group.
Divisor formulas, exact support definitions and the real plot all use the
original LMFDB coordinates, including the labelled points at infinity
:math:`\infty_{Y=\eta}` with :math:`Y=\lim y/x^3`.
The LMFDB link opens the database entry.

Choose a divisor in **Divisor support** to show its signed support on the
real locus. Filled markers have positive coefficients; open markers have
negative coefficients. Points at infinity are listed in the plot annotation.
Nonreal points retain their exact coordinates in **Exact support and
coefficients**. Algebraic points are labelled :math:`P_i`, with their
number field polynomial and embedding index given below the plot.

The cut panels use a separate square model :math:`v^2=F(u)`, displayed
above the tree. Point images under this recorded coordinate change determine
the logarithmic centers.

**The tree in** :math:`\mathbb P^1` shows the complete integration cut,
including logarithmic centers used by other basis divisors in the same
computation. **The lift to** :math:`C` shows its lift.

Select a logarithmic center in the tree panel and an edge to inspect the path integrals.
Changing the center preserves the tree view's pan and zoom.
Clicking a lifted path also selects its positive-side sheet.
For :math:`\omega_j=u^{j-1}du/v`, the path table gives

.. math::

   I_j(e)=\int_e\omega_j,\qquad
   L_{j,a}^{\pm}(e)=\int_e\log^{\pm}(u-a)\,\omega_j.

The regulator matrix pairs the two indicated divisor classes
:math:`B(D_i),B(D_j)` and :math:`\Delta C` with
:math:`\delta_{11},\delta_{12},\delta_{22}`. Its row weights are
:math:`(1,2\operatorname{Re}(u),|u|^2)` and its determinant is :math:`D`.
Each divisor column is the sum of its logarithmic columns with the exact
saved coefficients, including multiplicities and negative terms.

The comparison box displays

.. math::

   \pi c^+(C)c^-(C)L''(M_{\mathrm{std}}(C),1)\approx qD.

**Digits to display** selects 10–100 significant digits from the saved
values. The path table and matrix display 10 significant digits; downloads
retain the decimal strings and original certified interval endpoints.
Rational factors are numerical guesses, shown when their relative
agreement exceeds 95 digits. Certified enclosures retain the experiment's
analytic assumptions; agreement does not prove exact equality.

Cuts and sheets
---------------

A Euclidean minimum spanning tree joins the branch values and every recorded logarithmic
center. Circles mark branch values and diamonds mark
logarithmic centers. Coincident coordinates give one vertex. A ray from the
rightmost vertex extends the cut to infinity. Each logarithmic branch cut
is the tree path from its center to infinity, highlighted in its colour.
The sheet view uses the projection

.. math::

   (\operatorname{Re}u,\operatorname{Im}u,
   \operatorname{arsinh}(\operatorname{Re}(e^{-0.55i}v))).

Projection crossings need not be intersections on the curve. The meshes
are display approximations.

Recompute the catalog
---------------------

From the repository root:

.. code-block:: sh

   sage -python -m hyperell_regulator.experiments.beilinson_report \
       --input results/beilinson_results_certified.json \
       --bound 1000 --digits 100 --explorer-inputs docs/explorer_inputs.json
   sage -python -m hyperell_regulator.experiments.small_conductor_explorer --jobs 4
   sage -pip install '.[visualization]'
   sage -python -m hyperell_regulator.experiments.small_conductor_surfaces --jobs 4

The report export reads the completed run without repeating its integrations
or L-values. The explorer checks the saved divisor matrices using fresh
lower-precision integrals on the full shared cuts; its displayed numerical
values remain those from the successful run. The remaining commands install
Triangle and build display meshes. Source hashes ignore docstrings,
comments and formatting.

Sphinx copies ``docs/_extra/explorer`` into the published site.
Building the documentation does not repeat numerical integration.
KaTeX and fonts are served locally; their licenses are under
``explorer/vendor/``.

Explore another curve
---------------------

Compute a result, then select **Import local result**:

.. code-block:: sh

   sage -python -m hyperell_regulator.regulator.explorer \
       --F='x^5+(x+1)^2' --centers -1 0 --prec=30 --out=curve.json

The exporter accepts squarefree polynomials of degree 5 through 12 and
one through 128 logarithmic centers. Python/JSON requests also accept exact
algebraic centers as minimal-polynomial coefficients with a zero-based root
index, ordered by real then imaginary part. Approximate coordinates are
used only for display.
