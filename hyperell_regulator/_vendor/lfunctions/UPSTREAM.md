# Generic L-function backend

Source: https://github.com/edgarcosta/lfunctions

Pinned revision: `42014e1204f91e93f72b46945c4f5701b3ed7b1b`.
The upstream MIT license is retained in `LICENSE`.

This snapshot includes the values-only and Fourier error-bound corrections
documented in `docs/ball-arithmetic.rst`. In addition, the package uses the
entire sinc function at interpolation nodes and corrects the conditional
computation of the first-derivative interpolation error. It also implements
the previously unfinished conversion from the completed first derivative
to L'(s), using the logarithmic derivative of the gamma factor. The build never
fetches a moving upstream branch. No rank or zero search is required for the
certified value interface.
