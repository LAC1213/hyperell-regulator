Curve explorer
==============

`Open the explorer <explorer/index.html>`_ to inspect the Stokes computation
for the selected degree-zero divisors on a genus two curve
:math:`C:y^2+h(x)y=f(x)` in its original LMFDB coordinates. The catalog contains 41 successful comparisons from the 100-digit survey
below conductor 1000 and the individual curves ``3125.a.3125.1`` and
``4096.b.65536.1``. Six curves carry the rational endomorphism algebra
label ``CM``. The full survey lists four and six integral Galois-trace
generators for the two additional curves, respectively. The two divisors
selected for each regulator matrix are shown as :math:`D_1,D_2`.

The viewer
----------

Choose a curve to see its equation and the two divisors used for its
regulator matrix. The LMFDB link beside the selector opens the database
entry. The full survey lists packet torsion orders and rational
endomorphism algebras. Packet order is the exponent
:math:`N=\operatorname{lcm}_P\operatorname{ord}(P-P_0)` in the Jacobian, rather than
the number of points or the order of the ambient Jacobian torsion group.
Divisor formulas and exact support definitions use the original LMFDB
coordinates, including the labelled points at infinity
:math:`\infty_{\eta}` with :math:`\eta=\lim y/x^3`.
Quadratic coordinates are written with radicals. When a point is written
using an indexed root of a polynomial, the indices start at zero and
the roots are ordered by real part, then imaginary part.

The first panel displays the original LMFDB equation. The cut panels use
the completed square model :math:`v^2=f(u)`, displayed above the tree.
The transformed point images determine the logarithmic centers. This is
the paper's square model :math:`y^2=f(x)` with coordinates renamed to
distinguish it from the original LMFDB equation.

Integrality of the divisor classes
----------------------------------

Choose a bad reduction prime :math:`p` in **Integrality of α(Dᵢ)**.
The selector lists only the bad primes of the chosen curve and updates
the graph when the selection changes. The viewer uses exported exact
local models over a finite extension. When the support field has several
places above :math:`p`, the place selector chooses a local specialization.

The displayed graph is a unit-edge resolution of the verified metric
model after finite base change. Chains of rational components resolve
node thicknesses, and a further quadratic base change separates a
self-node when needed. Component indices from the marked metric model
are retained. The exported reduction data record the local-model field,
valuation normalization and resolution scale.
The graph need not be the minimal model over :math:`\mathbb Q_p`.

The graph lists the selected divisor's points and signed coefficients
beside the components where they specialize. The displayed genus is the
genus of the component's normalization. The panel lists unresolved
specializations below the graph. Unsupported
wild graphs are labelled **Graph not implemented** in the selector and
panel. Computations exceeding the export budget are labelled separately
as **Graph computation unavailable**. Neither case displays a guessed graph.

Enter rational charges :math:`f_Z` in the component table and select
**Solve** to solve :math:`Lu+f=0`. Calculations use exact fractions
and fix :math:`u_{Z_0}=0`. A solution exists precisely when the total
charge is zero. **Use divisor charge** supplies :math:`f=\sum_P a_PZ_P`.
It is enabled only when
all required point specializations are known.

**Local obstruction** displays the tuple from the paper.

.. math::

   \operatorname{ob}_{\mathfrak p}(u)
   =\bigl(Q_u,(u_Z)_{g(Z)\geq 1}\bigr),\qquad
   Q_u(\gamma,\delta)=\sum_{e=(Z,W)}(u_Z+u_W)\gamma_e\delta_e.

The cycle matrix uses the cycle basis in the exported reduction data. The panel also
tests membership in :math:`\mathbb Q\operatorname{ob}_{\mathfrak p}(\mathbf1)`;
vanishing of the raw tuple is sufficient but is not required. This
membership is independent of adding a constant to the potential. It is
a local condition for the entered charge, not a global integrality
certificate. Select a graph vertex to focus its charge input.

Regulator cut and matrix
------------------------

**The tree** :math:`\Gamma\subset\mathbb P^1` shows the complete integration cut,
including logarithmic centers used by other basis divisors in the same
computation. **The lift** :math:`T\subset C` shows its lift.

Select a logarithmic center in the tree panel and an edge to inspect the path integrals.
The **Path integrals** panel can be collapsed. Selecting an edge reopens it.
Changing the center preserves the tree view's pan and zoom.
Clicking a lifted path also selects its positive-side sheet.
For :math:`\omega_k=u^{k-1}du/v`, the path table uses the paper's notation

.. math::

   I_{n,j}^{(k)}=\int_{e_{n,j}^+}\omega_k,\qquad
   L_{n,j}^{(k)}=\int_{e_{n,j}^+}\log^+(u-a)\,\omega_k.

Here :math:`n` labels the positive-side sheet and :math:`j` the tree edge.
The gluing is :math:`e_{n,j}^-=e_{\tau_j(n),j}^+`.
The viewer stores edges oriented toward infinity; the orientation selector
can reverse them, negating both integrals. The negative logarithm-side
option displays the integral with :math:`\log^-` explicitly, and the jump
is displayed as :math:`\log^--\log^+` for the stored oriented cut.

The indicated divisors define the integral classes
:math:`\alpha(D_i),\alpha(D_j)`. The displayed numerical matrix is the
paper's plane-integral matrix :math:`R(D_i,D_j)`, rather than the matrix
of Poincaré pairings before normalization. Write

.. math::

   G_D(z)=\sum_{P\text{ finite}}a_P\log|z-u(P)|,
   \qquad D=\sum_P a_P P.

Then, in the integration chart,

.. math::

   R(D_i,D_j)=\int_{\mathbb C}\frac{1}{|f(z)|}
   \begin{pmatrix}1\\2\operatorname{Re}(z)\\|z|^2\end{pmatrix}
   \begin{pmatrix}G_{D_i}(z)&G_{D_j}(z)&1\end{pmatrix}
   \,\mathrm{d}x\,\mathrm{d}y,\qquad z=x+iy.

The rows follow :math:`\delta_{1,1},\delta_{1,2},\delta_{2,2}`;
the final column accounts for the Tate class :math:`\Delta C`.
The determinant is denoted :math:`\det R(D_i,D_j)`.
Each divisor column is the sum of its logarithmic columns with the exact
saved coefficients, including multiplicities and negative terms.

The comparison box displays

.. math::

   \pi c^+(C)c^-(C)L''(M_{\mathrm{std}}(C),1)
   \approx q\det R(D_i,D_j).

As in the paper's computational convention, :math:`c^+(C),c^-(C)`
are the signed maximal minors of :math:`\operatorname{Re}H` and
:math:`\operatorname{Im}H`, respectively, for the full-cycle period matrix
:math:`H` in the integration chart. The symbol :math:`D` always denotes a
divisor. The displayed rank concerns the examined packet support;
a rational-packet rank can be a lower bound for the full geometric packet.

**Digits to display** selects 10–100 significant digits from the saved
values. The path table and matrix display 10 significant digits; downloads
retain the decimal strings and original certified interval endpoints.
Values whose magnitude is less than half a unit in the selected decimal
precision display as zero, suppressing tiny numerical residuals.
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
       --bound 1000 --digits 100 --explorer-inputs docs/explorer_inputs.json \
       --extra-results docs/explorer_extra_results.json
   sage -python -m hyperell_regulator.experiments.small_conductor_explorer --jobs 4
   sage -pip install '.[visualization]'
   sage -python -m hyperell_regulator.experiments.small_conductor_surfaces --jobs 4
   sage -python -m hyperell_regulator.experiments.fibre_explorer --jobs 3 --seconds 45

The report export reads the completed run without repeating its integrations
or L-values. The explorer checks the saved divisor matrices using fresh
lower-precision integrals on the full shared cuts; its displayed numerical
values remain those from the successful run. The remaining commands install
Triangle and build display meshes. Source hashes ignore docstrings,
comments and formatting.
The fibre exporter separately computes the marked local models; it does
not repeat regulator or :math:`L`-value integrations. Its cached records
are keyed by the exact input support and model implementation.

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
