Interactive curve explorer
===========================

`Open the explorer <explorer/index.html>`_. Each curve has one branch diagram
with two sheets, shared by all its logarithmic weights. The catalog retains
the successful genus-two survey curves of conductor below 1300, excluding
conductor 1116, and their 100-digit reference determinants and L-values.

Using the viewer
----------------

Choose a curve and a logarithmic weight. Changing the weight updates the
path integrals and pairings without changing the geometry. Click a colored
path to select an edge and a positive-side sheet. The selectors also let
you choose the logarithm bank and reverse the path orientation. Drag the
plot to rotate; use its toolbar to zoom or pan. The sheet selector isolates
one sheet and both banks of its boundary. The mesh checkbox hides the
translucent surfaces while retaining paths.

For :math:`\omega_j=x^{j-1}dx/y`, the table displays

.. math::

   I_j(e)=\int_e\omega_j,\qquad
   L_{j,a}^{\pm}(e)=\int_e\log^{\pm}(x-a)\,\omega_j.

The cut includes all branch values, all selected logarithmic centers and
a ray to infinity. Logarithmic jumps can vary between edges:

.. math::

   L_{j,a}^{+}(e)-L_{j,a}^{-}(e)=\Delta_a(e) I_j(e),
   \qquad \Delta_a(e)\in 2\pi i\mathbb Z.

Both banks use the same geometric orientation toward infinity. To traverse
a sheet boundary, reverse the negative bank. Numerical values are displayed
to 10 significant digits; downloads retain the full exported decimal strings.

Geometry and pairings
---------------------

The single plot uses the projection

.. math::

   (\operatorname{Re}x,\operatorname{Im}x,
   \operatorname{arsinh}(\operatorname{Re}(e^{-0.55i}y))).

Projection crossings need not be curve intersections. Meshes use continuous
square-root factors on a triangulated slit domain. Their labels are checked
around mesh cycles and against both banks of the integration cut. Mesh
vertices are inset slightly from cuts. Meshes are display approximations
and do not enter numerical quadrature.

The plotted ray is truncated; its integrals include the entire infinite
tail. Finite endpoints are included in quadrature but approached by the
display samples. Float64 determines geometry and topology; guarded
arbitrary precision determines the integrals.

For the selected center a, the pairing table gives the plane regulator
:math:`R_{ij}`, its associated pairing :math:`-4R_{ij}/\pi`, and the
unweighted period vector. Genus-two vectors use weights
:math:`(1,\operatorname{Re}(x),|x|^2)`; double the middle entry for the
Beilinson determinant. For a finite common pole b, subtract its logarithmic
column from each marked-point column.

The survey table uses :math:`q=\pi c^+c^-L''(\mathrm{std},1)/D` and full
cycle periods. Rational guesses are shown when their saved agreement
exceeds 95 digits. Orders and component verdicts concern the algebraic
cycles; the numerical cover always has degree two. Precision targets and
convergence checks are not certified interval error bounds.

Recompute the catalog
---------------------

From the repository root:

.. code-block:: sh

   sage -python -m hyperell_regulator.experiments.beilinson_report --explorer-inputs docs/explorer_inputs.json
   sage -python -m hyperell_regulator.experiments.small_conductor_explorer
   sage -pip install '.[visualization]'
   sage -python -m hyperell_regulator.experiments.small_conductor_surfaces --jobs 4

The first command exports integration models, marked points and survey
references without constructing torsion functions. The second integrates
one shared cut per curve and verifies its determinant against the saved
reference. The third installs Triangle, the optional constrained
triangulation dependency. Mesh generation does not recompute integrals.

Sphinx copies the static files in ``docs/_extra/explorer`` to the published
``explorer/`` directory. Building the documentation does not run numerical
integration.

Explore another curve
---------------------

Compute a local result and select **Import local result**:

.. code-block:: sh

   sage -python -m hyperell_regulator.regulator.explorer \
     --F='x^5+(x+1)^2' --centers -1 0 --prec=30 --out=curve.json

The exporter accepts squarefree polynomials of degree 5 through 12 and
one through eight rational logarithmic centers. New exports use schema
version 2; old function-cover exports must be regenerated.

Preview documentation locally:

.. code-block:: sh

   sage -python -m sphinx -W --keep-going -b html docs docs/_build/html
   python3 -m http.server --directory docs/_build/html 8000
