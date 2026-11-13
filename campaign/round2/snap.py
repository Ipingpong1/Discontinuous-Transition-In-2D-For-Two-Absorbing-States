"""dnsim v1.5 snapshot records (--snap): u64 t, then Lx/64 x Ly u64 words, bit b of word x in row y =
cell (64 x + b, y).  Readers for the roughness measurements."""
import numpy as np


def lattice_size(prefix):
    for line in open(prefix + ".txt"):
        if line.startswith("# Lx"):
            w = line.split()
            return int(w[2]), int(w[4])
    raise ValueError(prefix)


def read_snaps(prefix, tmin=0):
    """Yields (t, grid) with grid a (Ly, Lx) uint8 array."""
    Lx, Ly = lattice_size(prefix)
    W = Lx // 64
    rec = 8 + 8 * W * Ly
    with open(prefix + ".snap", "rb") as f:
        while True:
            head = f.read(8)
            if len(head) < 8:
                return
            t = int(np.frombuffer(head, dtype=np.int64)[0])
            body = f.read(8 * W * Ly)
            if len(body) < 8 * W * Ly:
                return
            if t < tmin:
                continue
            bits = np.unpackbits(np.frombuffer(body, dtype=np.uint8), bitorder="little")
            yield t, bits.reshape(Ly, Lx)


def iface_mask(g):
    """Cells with at least one disagreeing Moore neighbour (periodic)."""
    m = np.zeros(g.shape, dtype=bool)
    for dy in (-1, 0, 1):
        for dx in (-1, 0, 1):
            if dy or dx:
                m |= g != np.roll(np.roll(g, dy, axis=0), dx, axis=1)
    return m


def column_edges(m, yc, half):
    """Outermost interface rows per column within rows yc - half .. yc + half - 1 (periodic): returns
    (lo, hi) as float arrays of length Lx (nan where the column has no interface cell), in the
    coordinate y - yc."""
    Ly, Lx = m.shape
    rows = (np.arange(yc - half, yc + half)) % Ly
    sub = m[rows]                                  # (2 half, Lx), row index k <-> y = yc - half + k
    k = np.arange(2 * half)[:, None]
    any_ = sub.any(0)
    lo = np.where(any_, np.where(sub, k, 10 ** 9).min(0), np.nan) - half
    hi = np.where(any_, np.where(sub, k, -1).max(0), np.nan) - half
    return lo.astype(float), hi.astype(float)


def height_diff(h, rs):
    """C(r) = <(h(x + r) - h(x))^2>_x on the periodic column index, for the lags rs."""
    out = []
    for r in rs:
        d = np.roll(h, -r) - h
        out.append(np.nanmean(d * d))
    return np.array(out)
