Interactive cut-system explorer
===============================

`Open the explorer <explorer/index.html>`_. It contains the successful 100-digit Beilinson regulator computations for
genus 2 curves of conductor strictly below 1300, with their two marked torsion
functions.
Calculations run locally in SageMath before publication; the published
GitHub Pages site serves only static files.

Using the viewer
----------------

Choose a curve and either marked torsion function. Click an edge in the
:math:`w=\varphi` plane, then a lift in the three-dimensional projection.
The selectors provide the same choices without using the plot. Drag the
lifted plot to rotate, and use its toolbar to zoom or pan. Switching tools
preserves the current view; only the **Home** button resets the camera.
The **Display sheets** selector shows all sheets or isolates one surface
with both banks of its boundary. The **Surface mesh** checkbox hides or
shows the surface while retaining the lifted paths. Path selection, sheet
filters, and changes of curve or function preserve the camera and zoom.

The path table shows, for every holomorphic form
:math:`\omega_a=x^{a-1}dx/y`,

.. math::

   I_a(e)=\int_e\omega_a,
   \qquad L_a^\pm(e)=\int_e\log^\pm(\varphi)\,\omega_a.

All numerical values displayed by the viewer are rounded to **10 significant
digits**. Exact polynomial coefficients, point coordinates and torsion orders
remain exact. The catalog retains the full 100-digit survey reference values.
Path integrals are independently exported at 15 digits and checked against
the survey determinant; their JSON retains the integration precision.
Conductor 1116 is retained in the report JSON and excluded from its LaTeX
tables and this explorer.

Sheets, banks and orientation
-----------------------------

The path selector numbers lifts by their positive-bank sheet. A geometric
edge bounds a positive-bank sheet :math:`S_i` and a negative-bank sheet
:math:`S_j`; both are shown in the path panel. On the same oriented lift,

.. math::

   L_a^+=L_a^-+\Delta\log(\varphi)I_a,
   \qquad \Delta\log(\varphi)=\pm2\pi i.

The jump and both bank labels come from the integration program. Both banks
initially point toward infinity. For a sheet's oriented boundary, traverse
the positive bank toward infinity and the negative bank toward zero. The
orientation control negates both ordinary and logarithmic integrals.

The base plot retains the actual polygonal detours of the cut system.
Unramified corners are combined into one logical edge by summing their
integrals with the corresponding sheet labels. The final ray is truncated
only in the plot; its displayed integral includes the entire infinite tail.
Finite display samples approach the endpoints, while quadrature includes
the endpoints. Three-dimensional paths are the projection

.. math::

   (\operatorname{Re}w,\operatorname{Im}w,
     \operatorname{arsinh}(\operatorname{Re}(e^{-0.55i}q))),

with :math:`q=y` for functions of :math:`x` alone and :math:`q=x` otherwise.
Projection crossings need not be intersections on the curve. Surface meshes
cover a finite window of the slit :math:`w`-plane. Their sheet labels are
continued from the same initial fibre as the integrator, without crossing
the cut. Both boundary labels are independently checked against the saved
paths. Triangle vertices approach cuts and branch points by tiny insets to
avoid singular display coordinates; the mesh is an approximation of the
projected surface. It does not change any integral.

Pairings and the determinant
----------------------------

The pairing panel shows the package regulator :math:`R_{ab}`, its associated
pairing :math:`-4R_{ab}/\pi`, and the unweighted period :math:`P_{ab}`. These
values belong to the exact function in the expanded curve/function panel.

For each marked point :math:`T` of order :math:`N`, the exported function
satisfies :math:`(\varphi_T)=N(T-Q)`, where :math:`Q` is the common pole
on the displayed integration model; it can be finite or at infinity.
The two marked zeros are labelled :math:`P` and :math:`R`, so the function
selector lists :math:`\varphi_P` and :math:`\varphi_R` with their orders.

If :math:`(x-x(T))^N=c\operatorname{Norm}(\varphi_T)`, its normalized
regulator column is

.. math::

   \frac{w_{ab}}{N}\left(2R_{ab}(\varphi_T)+\log|c|P_{ab}\right),
   \qquad (w_{11},w_{12},w_{22})=(1,2,1).

For a finite pole, replace :math:`x-x(T)` by
:math:`(x-x(T))/(x-x(Q))` in the norm identity.

The catalog determinant uses the columns for :math:`P`, :math:`R`, and
the periods multiplied by the same row weights, in that order. Thus the
determinant :math:`D` uses :math:`(1,2\operatorname{Re}(z),|z|^2)`. It is
checked against the current 100-digit survey determinant. The catalog uses
full cycle periods and :math:`q=\pi c^+c^-L''(\mathrm{std},1)/D`.
A missing numerical comparison is displayed as a dash. A dagger marks an
undecided component verdict, without asserting integrality. L-values and
relative comparison errors are retained from the report. The report stores
the zeros as ``P,Q`` and the pole as ``R``; the viewer labels them ``P,R``
and ``Q`` respectively.
Rational guesses are displayed when agreement exceeds 95 digits; otherwise
the computed quotient is displayed as a decimal. The relative error refers
to the saved rational guess.

Recompute the complete catalog
------------------------------

From the repository root, with the package installed in SageMath:

.. code-block:: sh

   sage -python -m hyperell_regulator.experiments.beilinson_report --explorer-inputs docs/explorer_inputs.json
   sage -python -m hyperell_regulator.experiments.small_conductor_explorer --jobs 4

Then generate the display meshes (this step does not reintegrate any edge):

.. code-block:: sh

   sage -pip install '.[visualization]'
   sage -python -m hyperell_regulator.experiments.small_conductor_surfaces --jobs 4 --allow-missing

The optional visualization dependency is `Triangle
<https://rufat.be/triangle/>`_, used to enforce all cut segments in the
triangulation even when critical values cluster closely. The published
viewer needs only the exported JSON files and the bundled JavaScript.

The function :math:`\varphi_P` for ``1147.a.35557.2`` has nearly coincident
cut segments at floating point precision. Its paths and numerical integrals
are available; its surface mesh is omitted because the sheet labels could
not be validated. The viewer identifies this limitation when it is selected.
Omit ``--allow-missing`` to require every mesh to validate during export.

The inputs in ``docs/explorer_inputs.json`` are generated from
``reports/beilinson_100_digits.json``. They preserve the successful survey
models, marked points, common poles, exact functions, torsion orders and
full-precision reference values. Use a fresh output directory when these
inputs change.
The exporter constructs the exact functions by Riemann--Roch and integrates
each cover once, retaining both ordinary and logarithmic edge integrals.
It checks sheet closure, edge traces, agreement of the periods from the two
covers, display continuation residuals and the reference determinants.
Raw curve/function integrals are reused when the geometric inputs, precision
and integration source hashes match. Determinant columns are normalized
separately; their source hashes preserve the distinction from the original
edge computations.

The JSON files and the plotting library are served from
``docs/_extra/explorer``. Sphinx copies them to ``explorer/`` in its HTML
output, and the existing documentation workflow publishes that output to
GitHub Pages. Numerical regeneration is a separate local operation, so a
documentation build does not repeat the integrations.

Explore another curve
----------------------

Compute a result locally and select **Import local result** in the viewer:

.. code-block:: sh

   sage -python -m hyperell_regulator.regulator.explorer \
     --F='x^5+(x+1)^2' --A='x+1' --B='1' --prec=15 --out=cut-system.json

The function is :math:`\varphi=A(x)+B(x)y`, with rational functions
:math:`A,B` over :math:`\mathbb Q`. The exporter verifies the single-zero,
single-pole divisor hypothesis before integration. Current export limits
are model degree 5 through 12 and cover degree at most 24.
The numerical algorithm uses arbitrary precision for integration and
float64 for initial topology and display continuation. Its digit target
and convergence checks do not constitute certified interval error bounds.

To preview the published files locally, build the documentation and serve it:

.. code-block:: sh

   sage -python -m sphinx -W --keep-going -b html docs docs/_build/html
   python3 -m http.server --directory docs/_build/html 8000

Visit ``http://localhost:8000/explorer/``. Use an HTTP server; browsers may
block JSON loading when the HTML file is opened directly from disk.

The vendored Plotly.js 2.9.0 bundle is distributed under its MIT license;
see ``explorer/PLOTLY-LICENSE.txt``. No external script or font service is
needed to use the viewer.
