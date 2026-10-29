hyperell-regulator
==================

Numerical tests of Beilinson's conjecture for the standard motive of genus 2
curves over the rationals, using SageMath.

The package computes periods, regulator integrals, standard :math:`L`-values,
and reduction and component data. Regulators, periods and :math:`L`-values
support both float64 and arbitrary precision.

For a curve :math:`C`, the standard motive is the rank-5 summand

.. math::

   \mathrm{std} = \ker\!\left(\bigwedge^2 h^1(C)
                  \longrightarrow h^2(C)\right).

Start with the installation guide and the worked example, or browse the API
reference generated from the package's docstrings.

.. toctree::
   :maxdepth: 2
   :caption: User guide

   installation
   quickstart
   experiments
   explorer-guide

.. toctree::
   :maxdepth: 2
   :caption: Reference

   api/modules

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

`Open the interactive cut-system explorer <explorer/index.html>`_ for all
sixteen curves and thirty-two torsion functions in the paper.
