"""
Coolant channel geometry: passage dimensions along the contour, the channel
count rule, and reading a channel geometry file (from scripts/channel_solver.py,
or hand-edited) for the wall temperature solver.
"""
import math
import pathlib
import numpy as np


def channel_geometry(x, r, N, w, h, t_w):
    """
    Passage dimensions at every station. w and h may be scalars or (m,) arrays.
        A_c   flow area
        D_h   hydraulic diameter 4 A / perimeter (Eq. 6.18)
        b     circumferential span per channel, land centre to land centre, at
              the channel floor: 2 pi (r + t_w) / N  (Eq. 6.39 discussion)
        land  b - w
        ds    wall-parallel segment length
        r, N  hot-gas wall radius and channel count, kept for the fin factor
        s     arc length along the wall from station 0 (injector face)
    """
    w = np.broadcast_to(np.asarray(w, dtype=float), r.shape).copy()
    h = np.broadcast_to(np.asarray(h, dtype=float), r.shape).copy()
    A_c = w * h
    D_h = 4 * A_c / (2 * (w + h))
    b = 2 * math.pi * (r + t_w) / N
    land = b - w
    ds = np.sqrt(np.gradient(x)**2 + np.gradient(r)**2)
    if land.min() <= 0:
        raise ValueError(f"{N} channels of width {w.max()*1e3:.2f} mm don't fit: "
                         f"land goes to {land.min()*1e3:.2f} mm at the throat")
    s = np.concatenate(([0.0], np.cumsum(np.hypot(np.diff(x), np.diff(r)))))
    return dict(w=w, h=h, A_c=A_c, D_h=D_h, b=b, land=land, ds=ds, r=np.asarray(r, dtype=float), N=N, s=s)


def set_station_height(chan, i, h):
    """Change the height at station i only, keeping A_c and D_h consistent (Eq. 6.18)."""
    w = chan["w"][i]
    chan["h"][i] = h
    chan["A_c"][i] = w * h
    chan["D_h"][i] = 4 * chan["A_c"][i] / (2 * (w + h))


def channels_from_land(Rt, w, t_w, land_min):
    """
    Step 2, milled-channel equivalent of Eq. 6.30: the most channels of width w
    that fit around the throat (floor radius Rt + t_w) with land >= land_min.
    """
    return int(2 * math.pi * (Rt + t_w) // (w + land_min))


def read_channel_geometry(path, x, r, tol_mm=0.01):
    """
    Channel geometry file -> (N, w, h) at the stations x, in m.

    Format (what util/export.write_channel_geometry writes):
        # ... any comment lines ...
        # N_channels: 115
        x_mm,r_mm,w_mm,h_mm
        0.0000,56.6240,0.5000,4.0000
        ...
    w and h are interpolated onto x. The file's wall radius r_mm must match the
    current contour within tol_mm: the channels only mean something on the wall
    they were made for. If the inputs changed the contour (O/F, Pc, CR, L*...),
    re-run the channel solver or edit a fresh file.
    """
    path = pathlib.Path(path)
    lines = path.read_text(encoding="utf-8").splitlines()
    N = None
    for l in lines:
        if l.startswith("#") and "N_channels:" in l:
            N = int(l.split("N_channels:")[1].strip())
    if N is None:
        raise ValueError(f"{path}: no '# N_channels: <n>' line")
    i_head = next(i for i, l in enumerate(lines) if l.strip() and not l.startswith("#"))
    cols = [c.strip() for c in lines[i_head].split(",")]
    data = np.loadtxt(lines[i_head + 1:], delimiter=",", ndmin=2)
    col = {name: data[:, cols.index(name)] * 1e-3 for name in ("x_mm", "r_mm", "w_mm", "h_mm")}

    xf, rf = col["x_mm"], col["r_mm"]
    tol = tol_mm * 1e-3
    if abs(xf[0] - x[0]) > tol or abs(xf[-1] - x[-1]) > tol:
        raise ValueError(f"{path}: covers x = {xf[0]*1e3:.2f}-{xf[-1]*1e3:.2f} mm, but the current "
                         f"contour is {x[0]*1e3:.2f}-{x[-1]*1e3:.2f} mm. Made for a different contour.")
    dr = np.max(np.abs(np.interp(x, xf, rf) - r))
    if dr > tol:
        raise ValueError(f"{path}: wall radius differs from the current contour by up to "
                         f"{dr*1e3:.3f} mm (limit {tol_mm} mm). Made for a different contour.")
    return N, np.interp(x, xf, col["w_mm"]), np.interp(x, xf, col["h_mm"])
