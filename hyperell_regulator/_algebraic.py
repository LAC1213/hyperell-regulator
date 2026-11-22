"""Exact algebraic-number records with a fixed complex embedding."""


def encode_algebraic(value):
    from sage.all import QQ, QQbar
    value = QQbar(value)
    if value in QQ:
        return {'rational': str(QQ(value))}
    polynomial = value.minpoly()
    roots = sorted(polynomial.roots(QQbar, multiplicities=False),
                   key=lambda z: (z.real(), z.imag()))
    return dict(polynomial=list(map(str, polynomial)), root_index=roots.index(value))


def decode_algebraic(record):
    from sage.all import PolynomialRing, QQ, QQbar
    if not isinstance(record, dict):
        return QQbar(QQ(record))
    if 'rational' in record:
        return QQbar(QQ(record['rational']))
    polynomial = PolynomialRing(QQ, 'a')(record['polynomial'])
    roots = sorted(polynomial.roots(QQbar, multiplicities=False),
                   key=lambda z: (z.real(), z.imag()))
    return roots[record['root_index']]


def complex_coordinates(value):
    """Approximate coordinates for geometry only; exact values remain encoded."""
    z = complex(value)
    return [z.real, z.imag]
