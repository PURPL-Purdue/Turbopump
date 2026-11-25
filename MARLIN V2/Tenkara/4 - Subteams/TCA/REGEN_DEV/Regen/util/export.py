"""
Files for CAD and for other runs. No physics here; the run scripts call these.
Every file gets the run's provenance stamp (util/provenance.py).

    write_contour_dxf       output/contour.dxf: wall, channel floor, channel top for NX
    write_height_csv        output/channel_heights.csv + NX point files (radial offsets)
    write_channel_geometry  output/channel_geometry.csv: N, w(x), h(x) on the wall they
                            were sized for, the input format of the wall temp solver
"""
import math
import numpy as np

from engine.contour import wall_angle
from util import provenance as prov


def write_contour_dxf(path, d, t_w=None, x=None, r=None, h=None):
    """
    Half-section contour for NX, mm, in the XY plane (X along the engine
    axis, Y = radius, x = 0 at the injector face).

    Every offset is perpendicular to the wall, so the wall is t_w thick and the
    channel h tall everywhere, which is what the 1D solver assumes.

    Layers:
        WALL           gas-side wall, exact: 2 lines, 2 arcs, and the bell as a
                       degree-2 spline with control points N, Q, E (the exact
                       Bezier the solver uses)
        CHANNEL_FLOOR  wall offset by t_w. Lines stay lines and arcs stay
                       concentric arcs (radius R - t_w). The sharp
                       cylinder->cone corner leaves a gap on the outside, which
                       is filled with an arc of radius t_w about the corner, so
                       every floor point is exactly t_w from the wall. The bell
                       floor is a fitted spline: the offset of a parabola isn't a
                       parabola.
        CHANNEL_TOP    wall offset by t_w + h(x) at each station (x, r, h), the
                       same corner fill, then one fitted spline. h varies, so it
                       has no exact form.
    """
    import ezdxf
    from ezdxf.math import fit_points_to_cad_cv

    mm = 1e3
    doc = ezdxf.new("R2010")
    doc.header["$INSUNITS"] = 4                       # 4 = millimetres
    msp = doc.modelspace()
    xt, Rc, Rt, x0 = d["x_t"], d["Rc"], d["Rt"], d["x0"]
    th_c = math.radians(d["theta_conv_deg"])
    (xN, rN), (xQ, rQ), (xE, rE) = d["N"], d["Q"], d["E"]
    P = lambda px, py: (px * mm, py * mm)

    def fitted_spline(pts, t0, t1, a):
        """Cubic spline through pts (m), with the end tangent directions t0, t1."""
        bs = fit_points_to_cad_cv([(px * mm, py * mm, 0) for px, py in pts],
                                  tangents=[(t0[0], t0[1], 0), (t1[0], t1[1], 0)])
        msp.add_spline(dxfattribs=a).apply_construction_tool(bs)

    def corner_fan(off, n=12):
        """Points on the arc of radius off about the cylinder->cone corner."""
        return [(x0 - off * math.sin(p), Rc + off * math.cos(p))
                for p in np.linspace(0, -th_c, n)]

    def add_offset_wall(layer, off):
        doc.layers.add(layer)
        a = dict(layer=layer)
        s, c = math.sin(th_c), math.cos(th_c)          # cone normal is (s, c)
        msp.add_line(P(0, Rc + off), P(x0, Rc + off), dxfattribs=a)             # cylinder
        if off > 0:                                                              # corner fill
            msp.add_arc(P(x0, Rc), off * mm, 90 - d["theta_conv_deg"], 90, dxfattribs=a)
        msp.add_line(P(x0 + off * s, Rc + off * c),
                     P(d["x1"] + off * s, d["r1"] + off * c), dxfattribs=a)      # convergent cone
        # Both arc centres sit above the throat, on the coolant side, so the
        # offset arcs keep the centre and shrink by off. 270 deg = throat.
        msp.add_arc(P(xt, Rt + d["R_up"]), (d["R_up"] - off) * mm,
                    270 - d["theta_conv_deg"], 270, dxfattribs=a)               # upstream arc
        msp.add_arc(P(xt, Rt + d["R_dn"]), (d["R_dn"] - off) * mm,
                    270, 270 + d["theta_n_deg"], dxfattribs=a)                  # downstream arc
        if off == 0:
            msp.add_open_spline([P(xN, rN), P(xQ, rQ), P(xE, rE)],
                                degree=2, dxfattribs=a)                         # bell, exact
        else:
            t = np.linspace(0, 1, 200)
            bx = (1 - t)**2 * xN + 2 * (1 - t) * t * xQ + t**2 * xE
            br = (1 - t)**2 * rN + 2 * (1 - t) * t * rQ + t**2 * rE
            dx = 2 * (1 - t) * (xQ - xN) + 2 * t * (xE - xQ)
            dr = 2 * (1 - t) * (rQ - rN) + 2 * t * (rE - rQ)
            L = np.hypot(dx, dr)
            pts = list(zip(bx - off * dr / L, br + off * dx / L))
            th_n, th_e = math.radians(d["theta_n_deg"]), math.radians(d["theta_e_deg"])
            fitted_spline(pts, (math.cos(th_n), math.sin(th_n)),
                          (math.cos(th_e), math.sin(th_e)), a)                  # bell, fitted

    add_offset_wall("WALL", 0.0)
    if t_w is not None:
        add_offset_wall("CHANNEL_FLOOR", t_w)
    if h is not None:
        doc.layers.add("CHANNEL_TOP")
        off = t_w + np.asarray(h)
        phi = wall_angle(x, d)
        px, pr = x - off * np.sin(phi), r + off * np.cos(phi)
        before = x <= x0
        pts = (list(zip(px[before], pr[before]))
               + corner_fan(float(np.interp(x0, x, off)))
               + list(zip(px[~before], pr[~before])))
        pts = [p for i, p in enumerate(pts)                                     # drop repeats
               if i == 0 or math.dist(p, pts[i - 1]) > 1e-9]
        fitted_spline(pts, (1.0, 0.0), (math.cos(phi[-1]), math.sin(phi[-1])),
                      dict(layer="CHANNEL_TOP"))

    # What made this file: DXF custom header properties (metadata, no geometry,
    # so NX import is unaffected). Read with ezdxf: doc.header.custom_vars
    for k, v in prov.stamp().items():
        doc.header.custom_vars.append(f"regen_{k}", str(v))
    doc.header.custom_vars.append("regen_inputs", "output/inputs_used.yaml")

    path.parent.mkdir(exist_ok=True)
    doc.saveas(path)
    prov.record(path)
    print(f"  Contour DXF written to {path}")


def write_height_csv(path, x, r, chan, t_w):
    """
    Channel height contour for NX. h is measured radially off the channel
    floor (outside of the hot wall, r + t_w); r_floor / r_top are included in
    case the floor needs checking against the CAD model. All mm.
    """
    mm = 1e3
    path.parent.mkdir(exist_ok=True)
    r_floor = r + t_w
    data = np.column_stack([x, chan["h"], r_floor, r_floor + chan["h"]]) * mm
    np.savetxt(path, data, delimiter=",", fmt="%.4f",
               header=f"# {prov.stamp_line()}\n# inputs: output/inputs_used.yaml\n"
                      "x_mm,h_mm,r_floor_mm,r_top_mm", comments="")
    prov.record(path)
    print(f"  Channel heights written to {path}")

    # NX Spline > Points from File: "X Y Z" per line, space separated, no header.
    # X along the engine axis, radius on Z, Y = 0 (XZ plane, matching the NX contour).
    for name, rad in (("channel_floor", r_floor), ("channel_top", r_floor + chan["h"])):
        pts = np.column_stack([x, np.zeros_like(x), rad]) * mm
        dat = path.with_name(f"{name}.dat")
        np.savetxt(dat, pts, delimiter=" ", fmt="%.4f")       # no stamp inline: NX reads every line
        prov.record(dat)
        print(f"  NX points written to {dat}")


def write_channel_geometry(path, x, r, N, chan):
    """
    Channel geometry on the wall it was sized for: x, wall radius r, width w and
    height h at every station, mm, plus the channel count. Read back with
    jacket.geometry.read_channel_geometry (the wall temp solver's geometry_file).
    Copy it and edit w_mm / h_mm to try a modified channel (rounded or smoothed
    heights, printer limits) without re-sizing.
    """
    mm = 1e3
    path.parent.mkdir(exist_ok=True)
    data = np.column_stack([x, r, chan["w"], chan["h"]]) * mm
    np.savetxt(path, data, delimiter=",", fmt="%.6f",
               header=f"# {prov.stamp_line()}\n# inputs: output/inputs_used.yaml\n"
                      f"# N_channels: {N}\n"
                      "x_mm,r_mm,w_mm,h_mm", comments="")
    prov.record(path)
    print(f"  Channel geometry written to {path}")
