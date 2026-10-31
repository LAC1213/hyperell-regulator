Interactive cut-system explorer
===============================

`Open the explorer <explorer/index.html>`_. It contains sixteen genus 2
curves of small conductor and thirty-two marked torsion functions.
Calculations run locally in SageMath before publication; the published
GitHub Pages site serves only static files.

Using the viewer
----------------

Choose a curve and either marked torsion function. Click an edge in the
:math:`w=\varphi` plane, then a lift in the three-dimensional projection.
The selectors provide the same choices without using the plot. Drag the
lifted plot to rotate, and use its toolbar to zoom and reset the view.
The **Display sheets** selector shows all sheets or isolates one surface
with both banks of its boundary. The **Surface mesh** checkbox hides or
shows the surface while retaining the lifted paths. Path selection preserves
the camera and zoom.

The path table shows, for every holomorphic form
:math:`\omega_a=x^{a-1}dx/y`,

.. math::

   I_a(e)=\int_e\omega_a,
   \qquad L_a^\pm(e)=\int_e\log^\pm(\varphi)\,\omega_a.

All numerical values displayed by the viewer are rounded to **10 significant
digits**. Exact polynomial coefficients, point coordinates and torsion orders
remain exact. Exported JSON retains the extra digits used for verification.

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
at infinity. In degree six, :math:`Q` is the point where
:math:`y/x^3\to1`; in odd degree it is the unique point at infinity.
The two marked zeros are labelled :math:`P` and :math:`R`, so the function
selector lists :math:`\varphi_P` and :math:`\varphi_R` with their orders.

If :math:`(x-x(T))^N=c\operatorname{Norm}(\varphi_T)`, its normalized
regulator column is

.. math::

   \frac{w_{ab}}{N}\left(2R_{ab}(\varphi_T)+\log|c|P_{ab}\right),
   \qquad (w_{11},w_{12},w_{22})=(1,2,1).

The catalog determinant uses the columns for :math:`P`, :math:`R`, and
the periods multiplied by the same row weights, in that order. Thus the
determinant :math:`D` uses :math:`(1,2\operatorname{Re}(z),|z|^2)`. It is
checked against the rescaled stored reference determinant. The catalog uses
full cycle periods and :math:`q=\pi c^+c^-L''(\mathrm{std},1)/D`.
All sixteen quotients have numerical comparisons. The dagger for
``388.a.776.1`` marks its still-undecided component condition at 2; its
L-value now uses the supported split nodal-model Euler factor. L-values
are retained as reference data by the viewer and exporter.

Recompute the complete catalog
------------------------------

From the repository root, with the package installed in SageMath:

.. code-block:: sh

   sage -python -m hyperell_regulator.experiments.small_conductor_explorer --jobs 4

Then generate the display meshes (this step does not reintegrate any edge):

.. code-block:: sh

   sage -pip install '.[visualization]'
   sage -python -m hyperell_regulator.experiments.small_conductor_surfaces --jobs 4

The optional visualization dependency is `Triangle
<https://rufat.be/triangle/>`_, used to enforce all cut segments in the
triangulation even when critical values cluster closely. The published
viewer needs only the exported JSON files and the bundled JavaScript.

The inputs in ``reports/small_conductor_explorer_inputs.json`` preserve the
sixteen models, their marked points and torsion orders, and the reference
values used for numerical checks.
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
