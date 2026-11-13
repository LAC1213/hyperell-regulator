Interactive curve explorer
===========================

`Open the explorer <explorer/index.html>`_. Each curve has a two-dimensional spanning tree
and a three-dimensional sheet view, shared by all its logarithmic weights. The catalog retains
the successful genus-two survey curves of conductor below 1300, excluding
conductor 1116, and their 100-digit reference determinants and L-values.

Using the viewer
----------------

Choose a curve and a logarithmic weight. Changing the weight updates the
path integrals and pairings without changing the geometry. Click an edge
in the two-dimensional tree to select it, or a colored path in the sheet
view to select both an edge and a positive-side sheet. The selected tree
edge and both of its lifted paths are highlighted and follow the edge selector. The selectors also let
you choose the logarithm bank and reverse the path orientation. Drag the
plot to rotate; use its toolbar to zoom or pan. The sheet selector isolates
one sheet and both banks of its boundary. The mesh checkbox hides the
translucent surfaces while retaining paths. Switching between pan, zoom and
rotation, changing curves, or updating a selection preserves the camera.
Use **Home** to reset it; double-clicking does not reset the camera.

For :math:`\omega_j=x^{j-1}dx/y`, the table displays

.. math::

   I_j(e)=\int_e\omega_j,\qquad
   L_{j,a}^{\pm}(e)=\int_e\log^{\pm}(x-a)\,\omega_j.

The cut is a Euclidean minimum spanning tree of all branch values and
selected logarithmic centers, with an outward ray to infinity at the
rightmost vertex. Different branches share a single diagram. Logarithmic jumps can vary between edges:

.. math::

   L_{j,a}^{+}(e)-L_{j,a}^{-}(e)=\Delta_a(e) I_j(e),
   \qquad \Delta_a(e)\in 2\pi i\mathbb Z.

Both banks use the same geometric orientation toward infinity. To traverse
a sheet boundary, reverse the negative bank and follow the cyclic edge
order at each junction. Downloads include this boundary order as
``boundary_word`` entries ``[edge_index, orientation]``. These edge indices
and the vertex indices in ``edges[].vertices`` are zero based; a null vertex
is infinity. The viewer displays edge numbers starting at one. Numerical values are displayed
to 10 significant digits; downloads retain the full exported decimal strings.

Geometry and pairings
---------------------

The two-dimensional plot shows the MST in the x-plane with equal axis
scales. Circles mark branch values, diamonds mark logarithmic centers, and
the dashed ray continues to infinity. Edge numbers match the integral
selectors. The three-dimensional sheet view uses the projection

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

The table recomputes determinants and full periods with the current MST
code at a 100-digit target. L-values are reused from the cache after the
current precision and coefficient-count checks. Quotients and their
relative agreement are recomputed from these ingredients. The table uses
:math:`q=\pi c^+c^-L''(\mathrm{std},1)/D` and full
cycle periods. Rational guesses are shown when their recomputed agreement
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
one shared cut per curve, recomputes the 100-digit numerical table, and
verifies its determinants against the saved reference. Numerical source
hashes and L-value cache checks are included in the catalog. The third
installs Triangle, the optional constrained
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
