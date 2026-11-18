r"""Reduction and the integrality criterion for hyperelliptic curves over `\QQ`.

At odd primes we use cluster pictures as in Dokchitser--Dokchitser--
Maistret--Morgan, Theorems 1.8 and 8.5. At 2 we use Liu's reduction data.
Component tests refer to the stable model; linking chains are contracted
to nodes. Component tests at 2 are undecided.

EXAMPLES::

    sage: from hyperell_regulator.reduction import HyperellipticCurveWithReduction
    sage: R.<x> = ZZ[]
    sage: C = HyperellipticCurveWithReduction(f=x^5 + x^3 + 1)
    sage: C.bad_primes()
    [2, 53, 61]
    sage: C.is_semistable_at(53)
    True
    sage: [(c.kind, c.genus) for c in C.components(53)]
    [('principal', 1)]
    sage: C.same_component((0, 1), (0, -1), 53)["answer"]
    True

    sage: from contextlib import redirect_stdout
    sage: from io import StringIO
    sage: output = StringIO()
    sage: with redirect_stdout(output):
    ....:     C.report(P=(0, 1), Q=(0, -1))
    sage: 'conductor' in output.getvalue()
    True
"""

import re

from sage.all import (
    Graph,
    HyperellipticCurve,
    PolynomialRing,
    QQ,
    Qp,
    SageObject,
    ZZ,
    ascii_art,
    factor,
    gcd,
    genus2reduction,
    prod,
)
from sage.schemes.hyperelliptic_curves.constructor import (
    HyperellipticCurve as _sage_hyperelliptic_curve)
from sage.structure.dynamic_class import dynamic_class

from sage_cluster_pictures.cluster_pictures import Cluster


# ---------------------------------------------------------------------------
# Parsing Liu / PARI ``genus2red`` local data
# ---------------------------------------------------------------------------

_NU_TYPE_RE = re.compile(r"reduction at p:\s*(.*?)\s*page", re.S)
_F_RE = re.compile(r"\bf\s*=\s*(\d+)")
_STABLE_RE = re.compile(r"\(potential\) stable reduction:\s*\(([^)]*)\)")


def _parse_local_data(text):
    r"""Liu's stable type, Namikawa--Ueno type and conductor exponent.

    Missing values are ``None``.
    """
    out = {"stable_type": None, "nu_type": None, "conductor_exponent": None,
           "raw": text}
    m = _STABLE_RE.search(text)
    if m:
        out["stable_type"] = m.group(1).strip()
    m = _NU_TYPE_RE.search(text)
    if m:
        out["nu_type"] = m.group(1).strip()
    m = _F_RE.search(text)
    if m:
        out["conductor_exponent"] = ZZ(m.group(1))
    return out


# Liu's classification of the *potential* stable reduction of a genus 2 curve
# (see the documentation of ``ReductionData``), together with the toric rank of
# the semi-abelian reduction of the Jacobian that it forces:
#
#   (I)   smooth genus 2                                      -> t = 0
#   (II)  elliptic curve with one node                        -> t = 1
#   (III) projective line with two nodes                      -> t = 2
#   (IV)  two projective lines meeting at three points        -> t = 2
#   (V)   two elliptic curves meeting at one point            -> t = 0
#   (VI)  elliptic curve + nodal line meeting at one point    -> t = 1
#   (VII) as (VI) with both components singular               -> t = 2
#
_POTENTIAL_TORIC_RANK = {"I": 0, "II": 1, "III": 2, "IV": 2,
                         "V": 0, "VI": 1, "VII": 2}


def potential_toric_rank(stable_type):
    r"""The potential toric rank from Liu's type (I)--(VII).

    EXAMPLES::

        sage: from hyperell_regulator.reduction import potential_toric_rank
        sage: potential_toric_rank("II")
        1
        sage: potential_toric_rank("V")
        0
    """
    if stable_type is None:
        return None
    key = stable_type.split(",")[0].strip()
    if key not in _POTENTIAL_TORIC_RANK:
        return None
    return ZZ(_POTENTIAL_TORIC_RANK[key])


def nu_type_is_semistable(label):
    r"""Whether the Namikawa--Ueno type has all component multiplicities equal to one.

    EXAMPLES::

        sage: from hyperell_regulator.reduction import nu_type_is_semistable
        sage: nu_type_is_semistable("[I{1-0-0}]")
        True
        sage: nu_type_is_semistable("[I{1-1-0}]")
        True
        sage: nu_type_is_semistable("[V]")
        False
        sage: nu_type_is_semistable("[II-II-0]")
        False
        sage: nu_type_is_semistable("[VIII-1]")
        False
        sage: nu_type_is_semistable("[I*{1-0-0}]")
        False
    """
    if label is None:
        return None
    if "*" in label:
        return False
    for m in re.finditer(r"[IVXK]+", label):
        if m.group(0) != "I":
            return False
        # a digit in front of a numeral is a multiplicity, e.g. "2I{0}-0"
        if m.start() > 0 and label[m.start() - 1].isdigit():
            return False
    return True


# ---------------------------------------------------------------------------
# Components of the special fibre
# ---------------------------------------------------------------------------

class SpecialFibreComponent(SageObject):
    r"""An irreducible component of a semistable special fibre.

    ``kind`` is ``"principal"`` or ``"chain"``. The latter is contracted in
    the stable model. ``sign`` distinguishes the two übereven components.
    """

    def __init__(self, label, kind, cluster=None, sign=None, genus=None,
                 equation=None):
        self.label = label
        self.kind = kind
        self.cluster = cluster
        self.sign = sign
        self.genus = genus
        self.equation = equation
        self.multiplicity = ZZ(1)

    def _repr_(self):
        s = "%s: %s, genus %s" % (self.label, self.kind, self.genus)
        if self.equation is not None:
            s += ", %s" % (self.equation,)
        else:
            s += ", P^1"
        return s


class PointReduction(SageObject):
    r"""The reduction of a rational point on the stable special fibre.

    ``is_smooth`` excludes nodes. ``point`` is a residue-field pair,
    ``"infinity"`` or ``None``; ``reason`` describes nonsmooth reduction.
    """

    def __init__(self, prime, is_smooth, component=None, point=None,
                 reason=None, cluster=None, y_datum=None):
        self.prime = prime
        self.is_smooth = is_smooth
        self.component = component
        self.point = point
        self.reason = reason
        self.cluster = cluster
        self._y_datum = y_datum
        self._raw = None        # the (x, y) that produced this, in K(R)

    def _repr_(self):
        if self.is_smooth:
            return "smooth point %s on %s (p=%s)" % (
                self.point, self.component.label, self.prime)
        return "singular point of the special fibre at p=%s (%s)" % (
            self.prime, self.reason)


# ---------------------------------------------------------------------------
# The main class
# ---------------------------------------------------------------------------

def HyperellipticCurveWithReduction(f, h=0, prec=250):
    r"""The Sage curve `y^2+h(x)y=f(x)` over `\QQ`, with reduction methods.

    ``f`` and ``h`` are integral; ``prec`` is the `p`-adic working precision.

    EXAMPLES::

        sage: from hyperell_regulator.reduction import HyperellipticCurveWithReduction
        sage: R.<x> = ZZ[]
        sage: C = HyperellipticCurveWithReduction(x^5 + x^3 + 1)
        sage: C
        Hyperelliptic Curve over Rational Field defined by y^2 = x^5 + x^3 + 1
        sage: C.genus()
        2
        sage: C.bad_primes()
        [2, 53, 61]
        sage: C.is_tame(53), C.is_semistable_at(53)
        (True, True)

        sage: C.jacobian()
        Jacobian of Hyperelliptic Curve over Rational Field defined by y^2 = x^5 + x^3 + 1

        sage: HyperellipticCurveWithReduction(x^7 - x^5 + x + 1).genus()
        3
    """
    S = PolynomialRing(ZZ, 'x')
    try:
        f = S(f)
        h = S(h)
    except TypeError:
        raise TypeError("f and h must be polynomials with integer coefficients")
    F = 4 * f + h ** 2
    if F.degree() < 3:
        raise ValueError("4f + h^2 has degree %s; a hyperelliptic curve needs "
                         "degree at least 3" % F.degree())
    if F.discriminant() == 0:
        raise ValueError("4f + h^2 is not separable: the curve is singular")
    RQ = PolynomialRing(QQ, 'x')
    H = _sage_hyperelliptic_curve(RQ(f), RQ(h))
    H.__class__ = dynamic_class("HyperellipticCurveWithReduction",
                                (HyperellipticReduction,), type(H))
    H.f, H.h, H.F = f, h, F
    H.prec = ZZ(prec)
    H._cache = {}
    return H


#: Short alias for :func:`HyperellipticCurveWithReduction`, for interactive use.
ClusterCurve = HyperellipticCurveWithReduction


class HyperellipticReduction:
    r"""Reduction methods for Sage hyperelliptic curves.

    ``F=4f+h^2`` defines the square model.
    """

    # -- basic objects ------------------------------------------------------

    def curve(self):
        r"""The underlying curve, ``self``."""
        return self

    def square_model(self):
        r"""The model `Y^2=4f+h^2`, isomorphic over `\ZZ[1/2]`."""
        RQ = PolynomialRing(QQ, 'x')
        return HyperellipticCurve(RQ(self.F))

    def to_square_model(self, P):
        r"""Send an affine point to the square model by `Y=2y+h(x)`.

        EXAMPLES::

            sage: from hyperell_regulator.reduction import HyperellipticCurveWithReduction
            sage: R.<x> = ZZ[]
            sage: C = HyperellipticCurveWithReduction(x^5 - x, x)
            sage: C.to_square_model((1, 0))
            (1, 1)
        """
        x0, y0 = self._coords(P)
        return (x0, 2 * y0 + self.f.parent()(self.h)(x0))

    def _coords(self, P):
        r"""Validate an affine rational point and return its coordinates."""
        if isinstance(P, (tuple, list)) and len(P) == 2:
            x0, y0 = QQ(P[0]), QQ(P[1])
        else:
            # (X : Y : Z) on Sage's plane model, which is the degree deg(F)
            # homogenisation y^2 z^{d-2} + ... ; setting Z = 1 gives x = X/Z,
            # y = Y/Z.
            if len(tuple(P)) != 3:
                raise ValueError("a point must be given by 2 or 3 coordinates")
            if QQ(P[2]) == 0:
                raise NotImplementedError(
                    "points at infinity are not supported; translate the "
                    "model so that the points of interest are affine")
            x0, y0 = QQ(P[0]) / QQ(P[2]), QQ(P[1]) / QQ(P[2])
        RQ = PolynomialRing(QQ, 'x')
        if y0 ** 2 + RQ(self.h)(x0) * y0 != RQ(self.f)(x0):
            raise ValueError("the point (%s, %s) is not on %s" % (x0, y0, self))
        return (x0, y0)

    # -- global reduction data ---------------------------------------------

    def reduction_data(self):
        r"""Liu's reduction data from ``genus2reduction``."""
        if self.genus() != 2:
            raise ValueError(
                "genus2reduction (Liu's algorithm) applies only to genus 2; "
                "this curve has genus %s.  Use the cluster picture methods, "
                "which work at every odd prime of tame reduction."
                % self.genus())
        if "reduction_data" not in self._cache:
            RQ = PolynomialRing(QQ, 'x')
            self._cache["reduction_data"] = genus2reduction(RQ(self.h), RQ(self.f))
        return self._cache["reduction_data"]

    def conductor(self):
        r"""
        The conductor of the Jacobian of `C`.
        """
        return ZZ(self.reduction_data().conductor)

    def minimal_discriminant(self):
        r"""The minimal discriminant of `C`."""
        return ZZ(self.reduction_data().minimal_disc)

    def bad_primes(self):
        r"""The prime divisors of the minimal discriminant.

        EXAMPLES::

            sage: from hyperell_regulator.reduction import HyperellipticCurveWithReduction
            sage: R.<x> = ZZ[]
            sage: HyperellipticCurveWithReduction(x^5 + x^3 + 1).bad_primes()
            [2, 53, 61]
        """
        if "bad_primes" not in self._cache:
            primes = None
            if self.genus() == 2:
                try:
                    primes = [q for q, _ in
                              ZZ(self.minimal_discriminant()).factor()]
                except Exception:
                    primes = None
            if primes is None:
                # No minimal discriminant available: start from the primes
                # dividing the discriminant of the model, a superset of the
                # bad primes, and sieve the odd ones with the cluster picture.
                D = self.F.discriminant() * self.F.leading_coefficient()
                primes = []
                for q, _ in ZZ(2 * D).factor():
                    if q == 2:
                        primes.append(q)      # undecidable here; keep it
                        continue
                    try:
                        if not self.cluster_picture(q).has_good_reduction(
                                Qp(q, self.prec)):
                            primes.append(q)
                    except NotImplementedError:
                        primes.append(q)      # wild: cannot rule it out
            self._cache["bad_primes"] = sorted(primes)
        return list(self._cache["bad_primes"])

    def local_data(self, p):
        r"""Liu's data at ``p``: ``stable_type``, ``nu_type``, ``conductor_exponent``, ``raw``."""
        p = ZZ(p)
        ld = self.reduction_data().local_data
        if p not in ld:
            return {"stable_type": None, "nu_type": "[I{0-0-0}]",
                    "conductor_exponent": ZZ(0), "raw": "good reduction"}
        return _parse_local_data(ld[p])

    # -- cluster pictures (odd p) -------------------------------------------

    def cluster_picture(self, p):
        r"""The cluster picture of `Y^2=4f+h^2` over `\QQ_p`, for odd ``p``.

        EXAMPLES::

            sage: from hyperell_regulator.reduction import HyperellipticCurveWithReduction
            sage: R.<x> = ZZ[]
            sage: C = HyperellipticCurveWithReduction(x^5 + x^3 + 1)
            sage: C.cluster_picture(53)
            Cluster with 5 roots and 4 children
        """
        p = ZZ(p)
        if p == 2:
            raise ValueError("cluster pictures require odd residue "
                             "characteristic; p = 2 is handled by Liu's "
                             "algorithm, see local_data(2)")
        key = ("cluster", p)
        if key not in self._cache:
            K = Qp(p, self.prec)
            RK = PolynomialRing(K, 'x')
            self._cache[key] = Cluster.from_polynomial(RK(self.F))
        return self._cache[key]

    # -- tameness -----------------------------------------------------------

    def _require_cluster_picture(self, p):
        r"""The cluster picture at an odd prime; raise if unavailable."""
        p = ZZ(p)
        if p == 2:
            raise NotImplementedError(
                "cluster pictures need odd residue characteristic")
        try:
            return self.cluster_picture(p)
        except NotImplementedError:
            raise NotImplementedError(
                "the splitting field of 4f + h^2 over Q_%s could not be built; "
                "sage_cluster_pictures refuses it as possibly wildly ramified, "
                "and is_tame(%s) = %s" % (p, p, self.is_tame(p)))

    def is_tame(self, p):
        r"""Whether wild inertia acts trivially on `H^1(C)` at ``p``.

        Return ``True``, ``False`` or ``None`` if undecided.

        EXAMPLES::

            sage: from hyperell_regulator.reduction import HyperellipticCurveWithReduction
            sage: R.<x> = ZZ[]
            sage: C = HyperellipticCurveWithReduction(x^6 - 27)
            sage: C.is_tame(3), C.is_semistable_at(3)
            (True, False)

            sage: HyperellipticCurveWithReduction(x^5 - x, x).is_tame(2)
            True
        """
        p = ZZ(p)
        if not p.is_prime():
            raise ValueError("%s is not prime" % p)

        # The conductor exponent of J bounds things both ways, and for genus 2
        # Liu gives it -- but NOT always at 2.  When PARI's genus2red cannot
        # determine the reduction at 2 it reports only the prime-to-2
        # conductor, and does so without any caveat, so v_2 of it is 0 by
        # construction rather than because J has good reduction there.  The
        # detectable signal is that the local data at 2 carries no
        # Namikawa-Ueno type and no conductor exponent; in that case f_2 is
        # simply unknown.  (Seen on 388.a.776.1, where PARI reports conductor
        # 97 while the true conductor is 388 = 4*97.)
        if self.genus() == 2:
            try:
                fp = ZZ(self.conductor()).valuation(p)
                if p == 2 and self.local_data(2)["conductor_exponent"] is None:
                    fp = None
            except Exception:
                fp = None
            if fp is not None:
                if fp == 0:
                    # J has good reduction, so rho_l is unramified
                    # (Neron-Ogg-Shafarevich), hence tame
                    return True
                if fp > 2 * self.genus():
                    # tame would force f_p = 2g - dim V^I <= 2g
                    return False

        # semistable (in particular good) reduction implies tame, for every p
        if self.genus() == 2 or p != 2:
            try:
                if self.is_semistable_at(p) is True:
                    return True
            except Exception:
                pass

        if p == 2:
            return None       # no criterion below applies in residue char 2

        if p > 2 * self.genus() + 1:
            return True

        K = Qp(p, self.prec)
        if all(ZZ(fac.degree()) % p != 0
               for fac, _ in PolynomialRing(K, 'x')(self.F).factor()):
            return True

        try:
            return ZZ(self.cluster_picture(p).n_wild()) == 0
        except (NotImplementedError, ValueError, ArithmeticError):
            return None

    def semistabilising_index(self, p):
        r"""A tame ramification index giving semistable reduction at ``p``.

        Return 1 exactly when `C` is already semistable over `\QQ_p`.

        EXAMPLES::

            sage: from hyperell_regulator.reduction import HyperellipticCurveWithReduction
            sage: R.<x> = ZZ[]
            sage: HyperellipticCurveWithReduction(x^6 - 27).semistabilising_index(3)
            4
            sage: HyperellipticCurveWithReduction(x^5 + x^3 + 1).semistabilising_index(53)
            1
        """
        p = ZZ(p)
        R = self._require_cluster_picture(p)
        if R.is_semistable(Qp(p, self.prec)):
            return ZZ(1)
        e = ZZ(R.roots()[0].parent().absolute_e())
        principal = [s for s in R.all_descendants() if s.is_principal()]
        k = ZZ(1)
        while True:
            m = e * k
            if all(m * s.depth() in ZZ and m * s.nu() / 2 in ZZ
                   for s in principal):
                return m
            k += 1

    # -- semistability ------------------------------------------------------

    def semistability_at(self, p):
        r"""Semistability at ``p``, with the method and reason.

        The dictionary entry ``semistable`` is ``True``, ``False`` or ``None``.

        EXAMPLES::

            sage: from hyperell_regulator.reduction import HyperellipticCurveWithReduction
            sage: R.<x> = ZZ[]
            sage: HyperellipticCurveWithReduction(x^5 - x, x).semistability_at(65563)["semistable"]
            True
        """
        p = ZZ(p)
        if not p.is_prime():
            raise ValueError("%s is not prime" % p)
        if p not in self.bad_primes():
            return {"prime": p, "semistable": True, "method": "good reduction",
                    "reason": "p does not divide the minimal discriminant",
                    "good": True}

        if self.genus() == 2:
            liu_ss, liu_reason = self._liu_semistability(p)
        else:
            liu_ss, liu_reason = None, None

        if p != 2:
            try:
                R = self.cluster_picture(p)
            except NotImplementedError:
                # ``Cluster.from_polynomial`` refuses wildly ramified
                # splitting fields.  Wild ramification means p | e, hence
                # e >= p >= 3 > 2, so [DDMM, Thm 1.8](1) already fails and the
                # curve is not semistable -- but the package's test for
                # wildness is only a sufficient condition for it, so we defer
                # to Liu's data instead of concluding.
                return {"prime": p, "semistable": liu_ss,
                        "method": "Liu/Namikawa-Ueno (cluster pictures "
                                  "unavailable: possibly wild extension)",
                        "reason": liu_reason, "liu": liu_ss,
                        "liu_reason": liu_reason}
            K = Qp(p, self.prec)
            ss = bool(R.is_semistable(K))
            reason = self._semistability_reason(R, K, ss)
            rec = {"prime": p, "semistable": ss, "method": "cluster picture",
                   "reason": reason, "good": bool(R.has_good_reduction(K)),
                   "liu": liu_ss, "liu_reason": liu_reason}
            if liu_ss is not None and liu_ss != ss:
                rec["warning"] = ("cluster pictures say semistable = %s but "
                                  "Liu's data say %s (%s)"
                                  % (ss, liu_ss, liu_reason))
            return rec

        # ---- p = 2: cluster pictures do not apply ------------------------
        if self.genus() != 2:
            return {"prime": p, "semistable": None, "method": "none",
                    "reason": "cluster pictures need odd residue "
                              "characteristic, and Liu's algorithm applies "
                              "only to genus 2 (this curve has genus %s)"
                              % self.genus()}
        return {"prime": p, "semistable": liu_ss, "method": "Liu/Namikawa-Ueno",
                "reason": liu_reason, "liu": liu_ss, "liu_reason": liu_reason}

    def _liu_semistability(self, p):
        r"""Semistability from Liu's conductor and Namikawa--Ueno type.

        Return ``(answer, reason)``, with ``answer=None`` if undecided.
        """
        ld = self.local_data(p)
        t = potential_toric_rank(ld["stable_type"])
        fexp = ld["conductor_exponent"]
        nu = ld["nu_type"]

        if fexp is not None and t is not None and fexp != t:
            return (False,
                    "the potential stable reduction has type (%s), forcing "
                    "toric rank %s, while the conductor exponent is f_%s = %s; "
                    "a semistable abelian surface has f_p equal to its toric "
                    "rank (Grothendieck), so C is not semistable at %s"
                    % (ld["stable_type"], t, p, fexp, p))

        by_nu = nu_type_is_semistable(nu)
        if by_nu is None:
            extra = ""
            if fexp is not None and t is not None:
                extra = ("; the necessary condition f_%s = %s = toric rank "
                         "does hold" % (p, fexp))
            return (None,
                    "PARI's genus2red returned no Namikawa-Ueno type at %s "
                    "(its algorithm is incomplete at 2), and cluster pictures "
                    "need odd residue characteristic; potential stable "
                    "reduction type (%s)%s"
                    % (p, ld["stable_type"], extra))
        return (by_nu,
                "the special fibre of the minimal regular model has "
                "Namikawa-Ueno type %s, whose components %s multiplicity one"
                % (nu, "all have" if by_nu else "do not all have"))

    def _semistability_reason(self, R, K, ss):
        r"""The failed clause of the cluster semistability criterion, if any."""
        if ss:
            return ("all proper clusters are inertia stable, the splitting "
                    "field has e <= 2, and every principal cluster has "
                    "integral depth and even nu")
        bad = []
        for s in R.all_descendants():
            if s.is_principal():
                if s.depth() not in ZZ:
                    bad.append("principal cluster of size %s has depth %s "
                               "which is not an integer" % (s.size(), s.depth()))
                elif s.nu() / 2 not in ZZ:
                    bad.append("principal cluster of size %s has nu = %s "
                               "which is not even" % (s.size(), s.nu()))
        try:
            Kr = R.roots()[0].parent()
            e = Kr.absolute_e() / gcd(Kr.absolute_e(), K.absolute_e())
            if e > 2:
                bad.append("the splitting field of 4f + h^2 has ramification "
                           "degree %s > 2 over Q_p" % e)
            for s in R.all_descendants():
                if s.is_proper() and s.inertia() != s:
                    bad.append("the cluster of size %s is not inertia stable"
                               % s.size())
                    break
        except Exception as exc:
            bad.append("(while checking the splitting field: %s)" % exc)
        return "; ".join(bad) if bad else "not semistable"

    def is_semistable_at(self, p):
        r"""
        Whether `C` is semistable at ``p`` (``None`` if undecided).
        """
        return self.semistability_at(p)["semistable"]

    def is_semistable(self):
        r"""Whether `C` is semistable at every prime; ``None`` if undecided.

        EXAMPLES::

            sage: from hyperell_regulator.reduction import HyperellipticCurveWithReduction
            sage: R.<x> = ZZ[]
            sage: HyperellipticCurveWithReduction(x^5 - x, x).is_semistable()
            True
        """
        answers = [self.is_semistable_at(p) for p in self.bad_primes()]
        if any(a is False for a in answers):
            return False
        if any(a is None for a in answers):
            return None
        return True

    def semistability(self):
        r"""
        The dictionary ``p -> semistability_at(p)`` over all bad primes.
        """
        return {p: self.semistability_at(p) for p in self.bad_primes()}

    # -- components of the special fibre -----------------------------------

    def components(self, p, model="stable"):
        r"""The special-fibre components at an odd bad prime.

        ``model="stable"`` retains principal components; ``model="regular"``
        also includes the linking chains. Return :class:`SpecialFibreComponent` objects.

        EXAMPLES::

            sage: from hyperell_regulator.reduction import HyperellipticCurveWithReduction
            sage: R.<x> = ZZ[]
            sage: C = HyperellipticCurveWithReduction(x^5 + x^3 + 1)
            sage: [(c.kind, c.genus) for c in C.components(53)]
            [('principal', 1)]
        """
        p = ZZ(p)
        if p == 2:
            raise NotImplementedError(
                "the components at p = 2 are not computed: cluster pictures "
                "need odd residue characteristic.  Liu's Namikawa-Ueno type "
                "at 2 is %r" % (self.local_data(2)["nu_type"],))
        ss = self.is_semistable_at(p)
        if ss is not True:
            raise ValueError("C is not semistable at %s, so the special fibre "
                             "is not a stable curve; the Namikawa-Ueno type is "
                             "%r" % (p, self.local_data(p)["nu_type"]))
        R = self.cluster_picture(p)
        if model == "stable":
            return self._principal_components(R)
        if model == "regular":
            return self._regular_components(R)
        raise ValueError('model must be "stable" or "regular"')

    def _cluster_names(self, R):
        r"""Label the top cluster ``R`` and its descendants ``s1``, ``s2``, etc."""
        names = {}
        i = 0
        for s in R.all_descendants():
            if not s.is_proper():
                continue
            if s.is_top_cluster():
                names[s] = "R"
            else:
                i += 1
                names[s] = "s%s" % i
        return names

    def _principal_components(self, R):
        names = self._cluster_names(R)
        out = []
        for s in R.all_descendants():
            if not s.is_principal():
                continue
            name = names[s]
            poly = self._component_polynomial(s)
            if s.is_ubereven():
                # Gamma_s : Y^2 = theta^2 G(X)^2 splits into two copies of P^1
                th = self._theta_bar(s)
                G = prod((PolynomialRing(poly.base_ring(), 'X').gen()
                          - self._red_cluster(s, c)) for c in s.children()
                         if c.is_twin() and c.relative_depth() == QQ(1) / 2)
                for eps in (1, -1):
                    if th is None:
                        eqn = ("Y = %ssqrt(%s) * (%s), conjugate over the "
                               "residue field"
                               % ("" if eps > 0 else "-",
                                  s.theta_squared().unit_part().residue(), G))
                    else:
                        eqn = "Y = %s" % (eps * th * G)
                    out.append(SpecialFibreComponent(
                        "Gamma_%s^%s" % (name, "+" if eps > 0 else "-"),
                        "principal", cluster=s, sign=eps, genus=ZZ(0),
                        equation=eqn))
            else:
                out.append(SpecialFibreComponent(
                    "Gamma_%s" % name, "principal", cluster=s, sign=None,
                    genus=ZZ(s.genus()), equation="Y^2 = %s" % poly))
        return out

    def _regular_components(self, R):
        r"""The components of the minimal regular special fibre from the dual graph."""
        principal = {}
        for comp in self._principal_components(R):
            key = comp.cluster if comp.sign is None else (comp.cluster, comp.sign)
            principal[key] = comp
        out = []
        for v in R.dual_graph().vertices(sort=False):
            if v in principal:
                out.append(principal[v])
            else:
                out.append(SpecialFibreComponent(
                    self._chain_label(v, self._cluster_names(R)), "chain",
                    cluster=None, sign=None, genus=ZZ(0), equation=None))
        return out

    @staticmethod
    def _chain_label(v, names):
        r"""A label for a chain vertex of the dual graph."""
        parts = []
        for u in v[:-1]:
            if isinstance(u, Cluster):
                parts.append("Gamma_%s" % names.get(u, "?"))
            elif isinstance(u, tuple):
                parts.append("Gamma_%s^%s" % (names.get(u[0], "?"),
                                              "+" if u[1] > 0 else "-"))
            else:
                parts.append("+" if u > 0 else "-")
        return "L(%s)[%s]" % (", ".join(parts), v[-1])

    def special_fibre(self, p):
        r"""The stable and regular components, dual graph, node count and toric rank at ``p``.

        EXAMPLES::

            sage: from hyperell_regulator.reduction import HyperellipticCurveWithReduction
            sage: R.<x> = ZZ[]
            sage: d = HyperellipticCurveWithReduction(x^5 + x^3 + 1).special_fibre(53)
            sage: d["toric_rank"]
            1
        """
        p = ZZ(p)
        R = self.cluster_picture(p)
        out = {"prime": p,
               "cluster_picture": R,
               "cluster_picture_ascii": ascii_art(R),
               "stable_components": self.components(p, "stable"),
               "toric_rank": ZZ(R.potential_toric_rank()),
               "conductor_exponent": ZZ(R.conductor_exponent()),
               "namikawa_ueno_type": self.local_data(p)["nu_type"],
               "regular_components": None,
               "dual_graph": None,
               "nodes": None,
               "tamagawa_number": None,
               "note": None}
        try:
            G = R.dual_graph()
            out["dual_graph"] = G
            out["nodes"] = ZZ(len(G.edges(sort=False)))
            out["regular_components"] = self.components(p, "regular")
        except ArithmeticError as exc:
            out["note"] = ("the chains of the minimal regular model were not "
                           "determined: %s" % exc)
        try:
            out["tamagawa_number"] = ZZ(R.tamagawa_number())
        except Exception:
            pass
        return out

    # -- reduction of points -----------------------------------------------

    @staticmethod
    def _valuation(a):
        r"""The valuation of ``a``, normalized by `v(p)=1`."""
        try:
            return a.normalized_valuation()
        except AttributeError:
            return QQ(a.valuation()) / a.parent().absolute_e()

    def _locate(self, R, x0):
        r"""Locate ``x0`` in the cluster tree.

        Return ``("component", s)``, ``("link", s, c, m)`` or ``("outside", R)``.
        A link contracts to a node in the stable model.
        """
        v = self._valuation
        if v(x0 - R.center()) < R.depth():
            return ("outside", R)
        s = R
        while True:
            nxt, lev = None, None
            for c in s.children():
                if not c.is_proper():
                    continue
                m = v(x0 - c.center())
                if m > s.depth():
                    nxt, lev = c, m
                    break
            if nxt is None:
                return ("component", s)
            if lev >= nxt.depth():
                s = nxt
                continue
            return ("link", s, nxt, lev)

    # The reduction map onto an übereven component is computed here because
    # ``Cluster.red`` cannot construct its pair of projective lines when the
    # component polynomial is constant.  Other components use the patched
    # ``Cluster.red`` described in the module requirements.  The formulas are
    # derived from
    #
    #     F(x)/\pi^{\nu_s} = \bar\theta_s^2 \prod_{s' child}(X - red(s'))^{|s'|},
    #     X = (x - z_s)/\pi^{d_s},
    #
    # which follows from `\nu_s = v(\theta_s^2) + |s| d_s`; the component is
    # `\Gamma_s : Y^2 = \bar\theta_s^2 G_s(X)` with `G_s` the product of
    # `X - red(s')` over the odd children and `(X - red(s'))^2` over the twins
    # of relative depth 1/2, so `Y = (y/\pi^{\nu_s/2}) /
    # \prod_{\delta_{s'} > 1/2} (X - red(s'))^{\lfloor |s'|/2 \rfloor}`.

    @staticmethod
    def _pi_pow(Kr, t):
        r"""`\pi^t` with `v(p)=1`; the uniformizer exponent is ``t*e``."""
        n = QQ(t) * Kr.absolute_e()
        if n not in ZZ:
            raise ValueError("exponent %s is not in the value group of %s"
                             % (t, Kr))
        return Kr.uniformiser() ** ZZ(n)

    def _red_x(self, s, xK):
        r"""The residue of `(x-z_{\mathfrak s})/\pi^{d_{\mathfrak s}}`."""
        Kr = xK.parent()
        a = (xK - Kr(s.center())) / self._pi_pow(Kr, s.depth())
        if a.valuation() > 0:
            return Kr.residue_field()(0)
        return a.residue()

    def _red_cluster(self, s, c):
        r"""The reduction of a child cluster ``c`` on ``s``."""
        Kr = s.roots()[0].parent()
        return self._red_x(s, Kr(c.roots()[0]))

    def _component_polynomial(self, s):
        r"""The polynomial `\bar\theta_{\mathfrak s}^2G_{\mathfrak s}(X)` defining the component."""
        # Empty products in the dependency can leave a residue field scalar
        # rather than a constant polynomial (notably for übereven clusters).
        k = s.roots()[0].parent().residue_field()
        return PolynomialRing(k, 'xL')(s.component_polynomial())

    def _theta_bar(self, s):
        r"""`\bar\theta_{\mathfrak s}` in the residue field, or ``None`` if nonsquare."""
        Kr = s.roots()[0].parent()
        th2 = Kr.residue_field()(s.theta_squared().unit_part().residue())
        return th2.sqrt() if th2.is_square() else None

    def _branch_value(self, s, X, Y):
        r"""The branch invariant `Y/G(X)=\pm\bar\theta` for an übereven cluster."""
        G = prod(X - self._red_cluster(s, c) for c in s.children()
                 if c.is_twin() and c.relative_depth() == QQ(1) / 2)
        return Y / G

    def _reduce_affine(self, s, xK, yK):
        r"""The residue coordinates on the principal component ``s``.

        Raise ``ArithmeticError`` if the image does not lie on the component.
        """
        if s.is_ubereven():
            X = self._red_x(s, xK)
            Y = (yK / self._pi_pow(xK.parent(), s.nu() / 2)).residue()
            Y *= prod((X - self._red_cluster(s, c)) ** (-(c.size() // 2))
                      for c in s.children()
                      if c.relative_depth() > QQ(1) / 2)
            if Y**2 != self._component_polynomial(s)(X):
                raise ArithmeticError("reduced point is not on the component")
            return (X, Y)
        try:
            pt = s.red((xK, yK))
        except ValueError as exc:
            raise ArithmeticError(str(exc))
        return (pt[0], pt[1])

    def _reduce_at_infinity(self, R, xK, yK):
        r"""The branch value at infinity for a point outside the top cluster.

        Return ``None`` in odd degree, where infinity is unique.
        """
        N = self.F.degree()
        if N % 2:
            return None
        Kr = xK.parent()
        c = R.leading_coefficient()
        V = yK / (self._pi_pow(Kr, self._valuation(Kr(c)) / 2)
                  * (xK - Kr(R.center())) ** (N // 2))
        if V.valuation() != 0:
            raise ArithmeticError("V = %s is not a unit at infinity" % V)
        return V.residue()

    @staticmethod
    def _infinity_child(s):
        r"""The principal child receiving infinity for a nonprincipal top cluster.

        Return ``None`` unless the top has one singleton and one proper child.
        """
        if not (s.is_top_cluster() and not s.is_principal() and s.is_even()):
            return None
        ch = s.children()
        if len(ch) != 2:
            return None
        proper = [c for c in ch if c.is_proper()]
        if len(proper) != 1:
            return None
        c = proper[0]
        return c if (c.is_principal() and c.is_odd()) else None

    def reduce_point(self, p, P):
        r"""Reduce an affine rational point on the stable fibre at an odd bad prime.

        Return a :class:`PointReduction`.

        EXAMPLES::

            sage: from hyperell_regulator.reduction import HyperellipticCurveWithReduction
            sage: R.<x> = ZZ[]
            sage: C = HyperellipticCurveWithReduction(x^5 + x^3 + 1)
            sage: C.reduce_point(53, (0, 1)).is_smooth
            True
        """
        p = ZZ(p)
        if p == 2:
            raise NotImplementedError(
                "reduction of points at p = 2 is not implemented (cluster "
                "pictures need odd residue characteristic)")
        R0 = self._require_cluster_picture(p)
        semistable = bool(R0.is_semistable(Qp(p, self.prec)))
        x0, y0 = self.to_square_model(P)
        R = self.cluster_picture(p)
        Kr = R.roots()[0].parent()
        xK, yK = Kr(x0), Kr(y0)

        loc = self._locate(R, xK)
        names = self._cluster_names(R)
        raw = (xK, yK)
        comps = self._principal_components(R) if semistable else []

        def pick(s, sign=None):
            if not semistable:
                # Without semistability over Q_p there is no special fibre over
                # F_p to label; the components live over the residue field of
                # the semistabilising extension.
                return SpecialFibreComponent(
                    "Gamma_%s" % names[s], "principal", cluster=s, sign=sign,
                    genus=ZZ(s.genus()), equation=None)
            return [c for c in comps if c.cluster == s and c.sign == sign][0]

        if loc[0] == "link":
            _, s, c, m = loc
            return PointReduction(
                p, False,
                reason="the x-coordinate lies strictly between the discs of "
                       "the clusters %s (depth %s) and %s (depth %s), at "
                       "level %s: the point reduces onto the chain of P^1's "
                       "joining them, i.e. to a node of the stable model"
                       % (names[s], s.depth(), names[c], c.depth(), m))

        if loc[0] == "outside":
            if not R.is_principal():
                c = self._infinity_child(R)
                if c is not None:
                    comp = pick(c)
                    r = PointReduction(p, True, component=comp,
                                       point="infinity", cluster=c)
                    r._raw = raw
                    return r
                return PointReduction(
                    p, False,
                    reason="the x-coordinate lies outside the disc of the top "
                           "cluster, which is not principal: the region at "
                           "infinity is part of a linking chain, so the point "
                           "reduces to a node of the stable model")
            try:
                V = self._reduce_at_infinity(R, xK, yK)
            except ArithmeticError as exc:
                return PointReduction(
                    p, False,
                    reason="could not reduce to infinity on Gamma_R: %s" % exc)
            if R.is_ubereven() and semistable:
                # V = +- theta, the same branch invariant as in the affine case
                comp = self._ubereven_component(R, V, comps)
                r = PointReduction(p, True, component=comp,
                                   point="infinity", cluster=R, y_datum=V)
                r._raw = raw
                return r
            comp = pick(R)
            r = PointReduction(p, True, component=comp, point="infinity",
                               cluster=R)
            r._raw = raw
            return r

        s = loc[1]
        if not s.is_principal():
            c = self._infinity_child(s)
            if c is not None:
                # the complement of c's disc maps to the single point at
                # infinity of Gamma_c, not onto a chain
                comp = pick(c)
                r = PointReduction(p, True, component=comp, point="infinity",
                                   cluster=c)
                r._raw = raw
                return r
            kind = ("twin" if s.is_twin()
                    else "cotwin" if s.is_cotwin()
                    else "not proper" if not s.is_proper()
                    else "non-principal")
            return PointReduction(
                p, False,
                reason="the x-coordinate lands in the %s cluster %s of size "
                       "%s: the corresponding piece of the special fibre is a "
                       "chain of P^1's, contracted to a node in the stable "
                       "model"
                       % (kind, names[s], s.size()))

        if not semistable:
            # C is tame but not semistable over Q_p: it becomes semistable over
            # a tamely ramified L, and the stable model lives there.  The
            # cluster combinatorics, and hence which component a point lands
            # on, are unchanged; only the depths and nu scale by e(L).  The
            # x-coordinate of the reduction is still computable over Q_p, since
            # red_s divides by an element of valuation exactly d_s.  The
            # y-coordinate is not (it needs pi^(nu/2), which may not exist over
            # Q_p), but it is never needed on its own -- only branch
            # *comparisons* between two points, in which that factor cancels.
            X = self._red_x(s, xK)
            r = PointReduction(p, True, component=pick(s), point=(X, None),
                               cluster=s)
            r._raw = raw
            return r

        try:
            X, Y = self._reduce_affine(s, xK, yK)
        except ArithmeticError as exc:
            return PointReduction(p, False,
                                  reason="could not reduce onto Gamma_%s: %s"
                                         % (names[s], exc))

        if s.is_ubereven():
            b = self._branch_value(s, X, Y)
            comp = self._ubereven_component(s, b, comps)
            r = PointReduction(p, True, component=comp, point=(X, Y),
                               cluster=s, y_datum=b)
            r._raw = raw
            return r
        r = PointReduction(p, True, component=pick(s), point=(X, Y), cluster=s)
        r._raw = raw
        return r

    def _ubereven_component(self, s, branch, comps):
        r"""The übereven component with branch value ``branch``."""
        th = self._theta_bar(s)
        if th is None:
            raise ArithmeticError(
                "theta^2 = %s is not a square in the residue field, so "
                "Gamma_s^+ and Gamma_s^- are conjugate and carry no rational "
                "point, yet a point reduced onto them"
                % s.theta_squared().unit_part().residue())
        sign = 1 if branch == th else -1
        return [c for c in comps if c.cluster == s and c.sign == sign][0]

    def same_component(self, P, Q, p):
        r"""Whether `P` and `Q` reduce smoothly on the same stable component at ``p``.

        Return ``prime``, ``answer``, the reductions ``P`` and ``Q``, and ``reason``.

        EXAMPLES::

            sage: from hyperell_regulator.reduction import HyperellipticCurveWithReduction
            sage: R.<x> = ZZ[]
            sage: C = HyperellipticCurveWithReduction(x^5 + x^3 + 1)
            sage: C.same_component((0, 1), (0, -1), 53)["answer"]
            True
        """
        p = ZZ(p)
        rP = self.reduce_point(p, P)
        rQ = self.reduce_point(p, Q)
        if not rP.is_smooth:
            return {"prime": p, "answer": False, "P": rP, "Q": rQ,
                    "reason": "P does not reduce to a smooth point: %s"
                              % rP.reason}
        if not rQ.is_smooth:
            return {"prime": p, "answer": False, "P": rP, "Q": rQ,
                    "reason": "Q does not reduce to a smooth point: %s"
                              % rQ.reason}
        if rP.cluster != rQ.cluster:
            return {"prime": p, "answer": False, "P": rP, "Q": rQ,
                    "reason": "P reduces to %s and Q reduces to %s"
                              % (rP.component.label, rQ.component.label)}
        if rP.cluster.is_ubereven():
            # Same übereven cluster: Gamma_s = Gamma_s^+ u Gamma_s^-, the two
            # branches Y = +- theta*G(X).  Decide whether the branches agree.
            same = self._same_branch(rP.cluster, rP, rQ)
            lab = rP.component.label.rsplit("^", 1)[0]
            if same is None:
                return {"prime": p, "answer": None, "P": rP, "Q": rQ,
                        "reason": "both reduce onto the übereven pair %s^+, "
                                  "%s^-, but C is not semistable at %s and one "
                                  "point is at infinity while the other is "
                                  "affine, so the two branch invariants cannot "
                                  "be compared without the semistabilising "
                                  "extension" % (lab, lab, p)}
            return {"prime": p, "answer": same, "P": rP, "Q": rQ,
                    "reason": "both reduce onto the übereven cluster's pair of "
                              "components %s^+, %s^-, on the %s branch"
                              % (lab, lab, "same" if same else "opposite")}
        return {"prime": p, "answer": True, "P": rP, "Q": rQ,
                "reason": "both reduce to smooth points of %s"
                          % rP.component.label}

    def _same_branch(self, s, rP, rQ):
        r"""Whether two smooth reductions lie on the same übereven component.

        Compare branch invariants, cancelling common ramification factors.
        Return ``None`` for an unresolved mixed affine/infinite pair.
        """
        if rP._y_datum is not None and rQ._y_datum is not None:
            return bool(rP._y_datum == rQ._y_datum)
        xP, yP = rP._raw
        xQ, yQ = rQ._raw
        infP = (rP.point == "infinity")
        infQ = (rQ.point == "infinity")
        if infP != infQ:
            return None
        if infP:
            N = self.F.degree()
            if N % 2:
                return True             # a single point at infinity
            z = xP.parent()(s.center())
            t = (yP / yQ) * ((xQ - z) / (xP - z)) ** (N // 2)
            return bool(t.residue() == 1)
        twins = [c for c in s.children()
                 if c.is_twin() and c.relative_depth() == QQ(1) / 2]
        XP, XQ = self._red_x(s, xP), self._red_x(s, xQ)
        GP = prod(XP - self._red_cluster(s, c) for c in twins)
        GQ = prod(XQ - self._red_cluster(s, c) for c in twins)
        return bool((yP / yQ).residue() * GQ == GP)

    def same_component_everywhere(self, P, Q):
        r"""Component tests keyed by bad prime; ``answer=None`` where unavailable."""
        out = {}
        for p in self.bad_primes():
            if p == 2:
                out[p] = {"prime": p, "answer": None, "P": None, "Q": None,
                          "reason": "not computed: cluster pictures need odd "
                                    "residue characteristic"}
                continue
            try:
                self._require_cluster_picture(p)
            except NotImplementedError as exc:
                out[p] = {"prime": p, "answer": None, "P": None, "Q": None,
                          "reason": str(exc)}
                continue
            out[p] = self.same_component(P, Q, p)
        return out

    # -- reporting ---------------------------------------------------------

    def report(self, P=None, Q=None):
        r"""Print reduction data and, if given, the component test for `P` and `Q`.

        EXAMPLES::

            sage: from hyperell_regulator.reduction import HyperellipticCurveWithReduction
            sage: R.<x> = ZZ[]
            sage: HyperellipticCurveWithReduction(x^5 - x, x).report(P=(0, 0), Q=(1, 0))  # random
        """
        print(self)
        print("  4f + h^2 = %s" % self.F)
        try:
            print("  conductor            = %s" % factor(self.conductor()))
            print("  minimal discriminant = %s" % factor(self.minimal_discriminant()))
        except Exception as exc:
            print("  (Liu's algorithm failed: %s)" % exc)
        bad = self.bad_primes()
        print("  bad primes           = %s" % bad)

        ss = self.is_semistable()
        print("")
        print("Semistable at all primes: %s" % ss)
        for p in bad:
            d = self.semistability_at(p)
            print("  p = %-8s semistable: %-5s  [%s]" % (p, d["semistable"],
                                                         d["method"]))
            print("      %s" % d["reason"])

        if ss is False:
            print("")
            print("Not semistable everywhere: skipping the special fibres.")
        else:
            for p in bad:
                print("")
                print("-" * 70)
                print("Special fibre at p = %s" % p)
                if p == 2:
                    ld = self.local_data(2)
                    print("  not computed (cluster pictures need p odd).")
                    print("  Liu: potential stable reduction type (%s), "
                          "Namikawa-Ueno type %s"
                          % (ld["stable_type"], ld["nu_type"]))
                    continue
                if self.is_semistable_at(p) is not True:
                    print("  C is not semistable at %s." % p)
                    continue
                d = self.special_fibre(p)
                print("  cluster picture:      %s"
                      % d["cluster_picture_ascii"])
                print("  Namikawa-Ueno type:   %s (Liu)"
                      % d["namikawa_ueno_type"])
                print("  toric rank %s, conductor exponent %s, "
                      "Tamagawa number %s"
                      % (d["toric_rank"], d["conductor_exponent"],
                         d["tamagawa_number"]))
                print("  components of the stable model (%s):"
                      % len(d["stable_components"]))
                for c in d["stable_components"]:
                    print("    %s" % c)
                if d["regular_components"] is None:
                    print("  minimal regular model: %s" % d["note"])
                else:
                    print("  components of the minimal regular model "
                          "(%s components, %s nodes):"
                          % (len(d["regular_components"]), d["nodes"]))
                    for c in d["regular_components"]:
                        print("    %s" % c)

        if P is None or Q is None:
            return
        print("")
        print("=" * 70)
        print("Reduction of P = %s and Q = %s" % (self._coords(P),
                                                  self._coords(Q)))
        for p, rec in self.same_component_everywhere(P, Q).items():
            print("  p = %-8s same smooth component: %s" % (p, rec["answer"]))
            print("      %s" % rec["reason"])
            for nm in ("P", "Q"):
                if rec[nm] is not None:
                    print("      %s -> %s" % (nm, rec[nm]))


# ---------------------------------------------------------------------------
# convenience
# ---------------------------------------------------------------------------

def genus2_semistable_report(f, h=0, P=None, Q=None, prec=250):
    r"""Construct the curve, print its reduction report and return it.

    EXAMPLES::

        sage: from hyperell_regulator.reduction import genus2_semistable_report
        sage: R.<x> = ZZ[]
        sage: C = genus2_semistable_report(x^5 - x, x, P=(0, 0), Q=(1, 0))  # random
    """
    C = HyperellipticCurveWithReduction(f, h, prec=prec)
    C.report(P=P, Q=Q)
    return C


def _liu_component_group_order(raw):
    r"""The geometric Néron component-group order from Liu's data, or ``None``."""
    m = re.search(r"page[^,]*,\s*(.*?),\s*f\s*=", raw, re.S)
    if not m:
        return None
    seg = m.group(1)
    if "H" in seg:
        return ZZ(4)
    order, found = ZZ(1), False
    for mm in re.finditer(r"\((\d+)\)(?:\^(\d+))?", seg):
        found = True
        n = max(ZZ(1), ZZ(mm.group(1)))
        k = ZZ(mm.group(2)) if mm.group(2) else ZZ(1)
        order *= n ** k
    return order if found else None


def self_test(verbose=True):
    r"""Check reduction, conductors and component groups against Liu's data.

    EXAMPLES::

        sage: from hyperell_regulator.reduction import self_test
        sage: self_test(verbose=False)  # long time
        True
    """
    S = PolynomialRing(ZZ, 'x')
    x = S.gen()
    tests = [
        (x**5 + x**3 + 1, S(0)),
        (x**6 + x**2 + 2, S(0)),
        (x**5 - x, x),
        (x**5 + 2*x**4 + 4, x + 1),
        (x**5 + 8, S(0)),
        (x**6 - 16, S(0)),
        (x*(x - 1)*(x - 2)*(x - 3)*(x - 4), S(0)),
        ((x**2 - 9)*(x**2 - 25)*(x**2 - 49), S(0)),
        ((x**2 - 5**4)*(x + 1)*(x + 2)*(x + 3), S(0)),
        ((x**4 - 5**8)*(x**2 + 2*x + 1 - 5**2), S(0)),
        (2*x**6 - 3*x**5 - 3*x**4 - 2*x**3 + 5*x + 1, S(2)),
        # two triples far apart: two genus 1 components and a chain
        (x*(x - 3**4)*(x + 3**4)*(x - 1)*(x - 1 - 3**4)*(x - 1 + 3**4), S(0)),
        # a long chain
        ((x - 1)*(x - 2)*(x - 3)*(x - 5**2)*(x - 5**5)*(x + 5**5), S(0)),
        # cotwin top cluster (child of size 2g)
        ((x - 1)*x*(x - 5**2)*(x - 2*5**2)*(x - 3*5**2), S(0)),
        ((x - 1)*x*(x - 5**4)*(x - 2*5**4)*(x - 3*5**4), S(0)),
        # cotwin that is not the top cluster
        ((x - 1)*(x - 5**2)*x*(x - 5**6)*(x - 2*5**6)*(x - 3*5**6), S(0)),
        ((x - 1)*(x - 5**4)*x*(x - 5**6)*(x - 2*5**6)*(x - 3*5**6), S(0)),
    ]
    failures = []
    for (f, h) in tests:
        C = HyperellipticCurveWithReduction(f, h, prec=300)
        for p in C.bad_primes():
            if p == 2:
                continue
            rec = C.semistability_at(p)
            if rec.get("liu") is not None and rec["liu"] != rec["semistable"]:
                failures.append((f, h, p, "semistability disagrees with Liu"))
            if rec["semistable"] is not True:
                continue
            d = C.special_fibre(p)
            fl = C.local_data(p)["conductor_exponent"]
            if fl is not None and d["conductor_exponent"] != fl:
                failures.append((f, h, p, "conductor exponent %s vs Liu %s"
                                 % (d["conductor_exponent"], fl)))
            G = d["dual_graph"]
            if G is None:
                continue
            nv = len(G.vertices(sort=False))
            ne = len(G.edges(sort=False))
            b1 = ne - nv + G.connected_components_number()
            gsum = sum(c.genus for c in d["regular_components"])
            if gsum + b1 != 2:
                failures.append((f, h, p, "genus formula: %s + %s != 2"
                                 % (gsum, b1)))
            if b1 != d["conductor_exponent"]:
                failures.append((f, h, p, "b_1 = %s but f_p = %s"
                                 % (b1, d["conductor_exponent"])))
            order = _liu_component_group_order(C.local_data(p)["raw"])
            if order is not None:
                H = Graph(multiedges=True)
                H.add_vertices(G.vertices(sort=False))
                H.add_edges([e[:2] for e in G.edges(sort=False)
                             if e[0] != e[1]])
                ntrees = (ZZ(1) if len(H.vertices(sort=False)) == 1
                          else ZZ(H.spanning_trees_count()))
                if ntrees != order:
                    failures.append((f, h, p, "component group %s vs Liu %s"
                                     % (ntrees, order)))
            # points
            RQ = PolynomialRing(QQ, 'x')
            fQ, hQ = RQ(f), RQ(h)
            pts = []
            for x0 in [QQ(a) for a in range(-6, 7)]:
                disc = hQ(x0) ** 2 + 4 * fQ(x0)
                if disc >= 0 and QQ(disc).is_square():
                    r = QQ(disc).sqrt()
                    for y0 in {(-hQ(x0) + r) / 2, (-hQ(x0) - r) / 2}:
                        pts.append((x0, y0))
            pts = pts[:6]
            for A in pts:
                rA = C.reduce_point(p, A)
                if rA.is_smooth and \
                        C.same_component(A, A, p)["answer"] is not True:
                    failures.append((f, h, p, "same_component not reflexive "
                                              "at %s" % (A,)))
                for B in pts:
                    if (C.same_component(A, B, p)["answer"]
                            != C.same_component(B, A, p)["answer"]):
                        failures.append((f, h, p, "same_component asymmetric "
                                                  "at %s, %s" % (A, B)))
        if verbose:
            print("checked %s" % HyperellipticCurveWithReduction(f, h))
    if verbose:
        print("")
        print("%s failure(s)" % len(failures))
        for t in failures:
            print("   ", t)
    return not failures


