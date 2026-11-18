Curve explorer
==============

`Open the explorer <explorer/index.html>`_ to inspect the Stokes computation
for the classes :math:`\alpha_P,\alpha_Q` on a genus two curve
:math:`C:y^2=F(x)`. The catalog contains 21 successful survey curves of
conductor below 1300 with 100-digit reference data.

The viewer
----------

Choose a curve to see its equation, marked points :math:`P,Q`, torsion
orders and real locus. The LMFDB link opens its database entry.
**The tree in** :math:`\mathbb P^1` shows the cut on the base;
**The lift to** :math:`C` shows its lift.

Select a logarithmic center in the tree panel and an edge to inspect the path integrals.
Changing the center preserves the tree view's pan and zoom.
Clicking a lifted path also selects its positive-side sheet.
For :math:`\omega_j=x^{j-1}dx/y`, the path table gives

.. math::

   I_j(e)=\int_e\omega_j,\qquad
   L_{j,a}^{\pm}(e)=\int_e\log^{\pm}(x-a)\,\omega_j.

The regulator matrix pairs :math:`\alpha_P,\alpha_Q,\Delta C` with
:math:`\delta_{11},\delta_{12},\delta_{22}`. Its row weights are
:math:`(1,2\operatorname{Re}(x),|x|^2)` and its determinant is :math:`D`.
A finite common pole contributes its logarithmic column with a minus sign
to both marked-point columns.

The comparison box displays

.. math::

   \pi c^+(C)c^-(C)L''(M_{\mathrm{std}}(C),1)=qD.

**Digits to display** selects 10–100 significant digits from the saved
values. The path table and matrix display 10 significant digits; downloads
retain the full decimal strings. Rational factors are numerical guesses,
shown when their relative agreement exceeds 95 digits.

Cuts and sheets
---------------

A Euclidean minimum spanning tree joins the branch values and marked
:math:`x`-coordinates. Circles mark branch values and diamonds mark
logarithmic centers. Coincident coordinates give one vertex. A ray from the
rightmost vertex extends the cut to infinity. Each logarithmic branch cut
is the tree path from its center to infinity, highlighted in its colour.
The sheet view uses the projection

.. math::

   (\operatorname{Re}x,\operatorname{Im}x,
   \operatorname{arsinh}(\operatorname{Re}(e^{-0.55i}y))).

Projection crossings need not be intersections on the curve. The meshes
are display approximations.

Recompute the catalog
---------------------

From the repository root:

.. code-block:: sh

   sage -python -m hyperell_regulator.experiments.beilinson_report --explorer-inputs docs/explorer_inputs.json
   sage -python -m hyperell_regulator.experiments.small_conductor_explorer
   sage -pip install '.[visualization]'
   sage -python -m hyperell_regulator.experiments.small_conductor_surfaces --jobs 4

These commands export the survey models, integrate their common cuts,
install Triangle and build display meshes. Determinants are checked against
the survey references. Source hashes ignore docstrings, comments and formatting.

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
one through eight rational logarithmic centers.
