"""Small rational linear algebra, deliberately independent of numerical libraries."""

from fractions import Fraction as F


class Invalid(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise Invalid(message)


def integer(x):
    return type(x) is int  # JSON booleans are not integers in our schema.


def transpose(a):
    return [list(t) for t in zip(*a)]


def eye(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def mul(a, b):
    return [[sum(x * y for x, y in zip(row, col)) for col in zip(*b)] for row in a]


def mv(a, v):
    return [sum(x * y for x, y in zip(row, v)) for row in a]


def rref(a):
    a = [[F(x) for x in row] for row in a]
    if not a:
        return a, []
    pivots = []
    r = 0
    for j in range(len(a[0])):
        k = next((k for k in range(r, len(a)) if a[k][j]), None)
        if k is None:
            continue
        a[r], a[k] = a[k], a[r]
        q = a[r][j]
        a[r] = [x / q for x in a[r]]
        for k in range(len(a)):
            if k != r and a[k][j]:
                q = a[k][j]
                a[k] = [x - q * y for x, y in zip(a[k], a[r])]
        pivots.append(j)
        r += 1
        if r == len(a):
            break
    return a, pivots


def rank(a):
    return len(rref(a)[1])


def inverse(a):
    n = len(a)
    require(n > 0 and all(len(row) == n for row in a), "inverse: nonsquare")
    e = eye(n)
    b, _ = rref([list(row) + e[i] for i, row in enumerate(a)])
    require([row[:n] for row in b] == e, "singular matrix")
    return [row[n:] for row in b]


def det(a):
    n = len(a)
    require(all(len(row) == n for row in a), "determinant: nonsquare")
    b = [[F(x) for x in row] for row in a]
    ans = F(1)
    for j in range(n):
        k = next((k for k in range(j, n) if b[k][j]), None)
        if k is None:
            return F(0)
        if k != j:
            b[k], b[j] = b[j], b[k]
            ans = -ans
        q = b[j][j]
        ans *= q
        for k in range(j + 1, n):
            t = b[k][j] / q
            for l in range(j + 1, n):
                b[k][l] -= t * b[j][l]
    return ans


def outer(v):
    return [[x * y for y in v] for x in v]


def symmetric_coords(a):
    # Upper triangle without doubling; an injective linear coordinate map.
    return [a[i][j] for i in range(len(a)) for j in range(i, len(a))]


def span_action(vectors, g):
    """Compute T_g on a greedy basis, and verify EVERY coordinate equality."""
    rays = [symmetric_coords(outer(v)) for v in vectors]
    indices = rref(transpose(rays))[1]
    basis = transpose([rays[i] for i in indices])
    rows = rref(transpose(basis))[1]
    left_inverse = inverse([basis[i] for i in rows])
    images = []
    for k in indices:
        # This is explicitly g A g^t, not a parity shortcut.
        a = mul(mul(g, outer(vectors[k])), transpose(g))
        images.append(symmetric_coords(a))
    images = transpose(images)
    action = mul(left_inverse, [images[i] for i in rows])
    require(mul(basis, action) == images, "action does not preserve W")
    return indices, action, det(action)


def rational_string(x):
    x = F(x)
    return str(x.numerator) if x.denominator == 1 else f"{x.numerator}/{x.denominator}"


def rational_matrix(a):
    return [[rational_string(x) for x in row] for row in a]
