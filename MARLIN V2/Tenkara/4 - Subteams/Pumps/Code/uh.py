"""
Meridional cross-section of a centrifugal impeller (incompressible) -- Plotly version.

Axial direction = x, radius = y (d/2), rotation axis at y = 0. Drawn 1:1 to scale.

Required parameters
    d2, b2, zE, d1, g1, eps_ds (deg), R_ds
Optional
    eps_ts (deg)  hub outlet angle            default = eps_ds
    R_ts          hub arc radius              default: hub arc starts at x = 0
    eps_ek (deg)  leading-edge angle          default = 45

Features: hover read-out (x, r, d), zoom/pan with locked 1:1 aspect, toggleable
legend entries, dimension arrows, optional mirrored lower half, HTML/PNG export,
and an optional 3D surface of revolution.

Requires:  pip install plotly            (PNG export also needs: pip install kaleido)
"""
import warnings
import numpy as np
import plotly.graph_objects as go


# =============================================================== geometry
def _arc(p0, R, phi0, phi1, n=200):
    """Arc: point(phi) = c + R*(sin phi, -cos phi), tangent (cos phi, sin phi).
    phi = 0 heads in +x (axial), phi = 90 deg heads in +y (radial). Starts at p0."""
    c = np.array([p0[0] - R * np.sin(phi0), p0[1] + R * np.cos(phi0)])
    phi = np.linspace(phi0, phi1, n)
    return np.column_stack((c[0] + R * np.sin(phi), c[1] - R * np.cos(phi)))


def shroud_profile(d1, g1, d2, zE, eps_ds, R_ds):
    """Straight g1 -> arc R_ds -> blending arc (radius solved) so that
    d1, g1, R_ds, zE, d2 and eps_ds are all honoured. Falls back to one arc
    (adapting R_ds and g1, with a warning) if that is impossible."""
    r1, r2 = d1 / 2, d2 / 2
    H = r2 - r1
    beta = np.pi / 2 - np.radians(eps_ds)

    def R2_of(a):
        return (H - R_ds * (1 - np.cos(a))) / (np.cos(a) - np.cos(beta))

    def res(a):
        return g1 + R_ds * np.sin(a) + R2_of(a) * (np.sin(beta) - np.sin(a)) - zE

    alphas = np.linspace(0.0, beta * (1 - 1e-4), 4000)
    with np.errstate(all="ignore"):
        R2, f = R2_of(alphas), res(alphas)
    ok = np.isfinite(R2) & (R2 > 0) & np.isfinite(f)

    alpha = None
    for i in range(len(alphas) - 1):
        if ok[i] and ok[i + 1] and f[i] * f[i + 1] <= 0:
            lo, hi = alphas[i], alphas[i + 1]
            for _ in range(80):
                mid = 0.5 * (lo + hi)
                lo, hi = (lo, mid) if res(lo) * res(mid) <= 0 else (mid, hi)
            alpha = 0.5 * (lo + hi)
            break

    if alpha is not None:
        a1 = _arc((g1, r1), R_ds, 0.0, alpha)
        a2 = _arc(a1[-1], R2_of(alpha), alpha, beta)
        return np.vstack([[0.0, r1], [g1, r1], a1, a2])

    R = H / (1 - np.cos(beta))
    g1_eff = zE - R * np.sin(beta)
    if g1_eff < 0:
        raise ValueError("Shroud parameters are inconsistent (zE too small for d1/d2/eps_ds).")
    warnings.warn(f"R_ds={R_ds}, g1={g1} incompatible with zE/d1/d2/eps_ds; "
                  f"using R_ds={R:.4g}, g1={g1_eff:.4g}.")
    return np.vstack([[0.0, r1], [g1_eff, r1], _arc((g1_eff, r1), R, 0.0, beta)])


def hub_profile(zE, b2, d2, eps_ts, R_ts):
    r2, xe = d2 / 2, zE + b2
    beta = np.pi / 2 - np.radians(eps_ts)
    R = xe / np.sin(beta) if R_ts is None else R_ts
    xs = xe - R * np.sin(beta)
    rn = r2 - R * (1 - np.cos(beta))
    if xs < -1e-9 or rn <= 0:
        raise ValueError("Hub arc does not fit; change R_ts / eps_ts.")
    arc = _arc((xs, rn), R, 0.0, beta)
    return np.vstack([[0.0, rn], arc]) if xs > 1e-9 else arc


def leading_edge(shroud, hub, eps_ek):
    y0 = shroud[0, 1]
    diff = hub[:, 1] - (y0 - np.tan(np.radians(eps_ek)) * hub[:, 0])
    idx = np.where(diff[:-1] * diff[1:] <= 0)[0]
    if len(idx) == 0:
        return None
    i = idx[0]
    t = diff[i] / (diff[i] - diff[i + 1]) if diff[i] != diff[i + 1] else 0.0
    return np.array([[0.0, y0], hub[i] + t * (hub[i + 1] - hub[i])])


# =============================================================== plotting
def _line_trace(xy, name, color="black", width=3, dash=None, showlegend=True, opacity=1.0):
    d = 2 * xy[:, 1]
    return go.Scatter(
        x=xy[:, 0], y=xy[:, 1], mode="lines", name=name, opacity=opacity,
        line=dict(color=color, width=width, dash=dash), showlegend=showlegend,
        customdata=d,
        hovertemplate=f"<b>{name}</b><br>x = %{{x:.2f}}<br>r = %{{y:.2f}}"
                      "<br>d = %{customdata:.2f}<extra></extra>",
    )


def _dim(fig, p0, p1, text, text_pos=None, color="#1f77b4"):
    """Double-headed dimension arrow from p0 to p1 with a label."""
    fig.add_annotation(x=p1[0], y=p1[1], ax=p0[0], ay=p0[1], axref="x", ayref="y",
                       xref="x", yref="y", showarrow=True, arrowhead=2, startarrowhead=2,
                       arrowwidth=1.2, arrowsize=1, arrowcolor=color, text="")
    mx, my = text_pos if text_pos else ((p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2)
    fig.add_annotation(x=mx, y=my, text=text, showarrow=False, font=dict(color=color, size=13),
                       bgcolor="rgba(255,255,255,0.8)")


def _extension(fig, p0, p1):
    fig.add_shape(type="line", x0=p0[0], y0=p0[1], x1=p1[0], y1=p1[1],
                  line=dict(color="gray", width=1, dash="dot"))


def build_geometry(d2, b2, zE, d1, g1, eps_ds, R_ds, eps_ts=None, R_ts=None, eps_ek=45.0):
    eps_ts = eps_ds if eps_ts is None else eps_ts
    shroud = shroud_profile(d1, g1, d2, zE, eps_ds, R_ds)
    hub = hub_profile(zE, b2, d2, eps_ts, R_ts)
    outlet = np.array([[zE, d2 / 2], [zE + b2, d2 / 2]])
    return shroud, hub, outlet, leading_edge(shroud, hub, eps_ek)


def plot_impeller(d2, b2, zE, d1, g1, eps_ds, R_ds,
                  eps_ts=None, R_ts=None, eps_ek=45.0,
                  dimensions=True, mirror=False, template="plotly_white"):
    """2D meridional section. Returns a plotly Figure."""
    shroud, hub, outlet, le = build_geometry(d2, b2, zE, d1, g1, eps_ds, R_ds,
                                             eps_ts, R_ts, eps_ek)
    r1, r2, xe = d1 / 2, d2 / 2, zE + b2
    fig = go.Figure()

    if mirror:  # faint lower half
        for xy in (shroud, hub, outlet):
            fig.add_trace(go.Scatter(x=xy[:, 0], y=-xy[:, 1], mode="lines", hoverinfo="skip",
                                     line=dict(color="lightgray", width=2), showlegend=False))

    fig.add_trace(_line_trace(shroud, "Shroud (DS)", "#d62728"))
    fig.add_trace(_line_trace(hub, "Hub (TS)", "#1f77b4"))
    fig.add_trace(_line_trace(outlet, "Outlet (b2)", "black"))
    if le is not None:
        fig.add_trace(_line_trace(le, "Leading edge", "#2ca02c", width=2))
    fig.add_trace(go.Scatter(x=[-0.2 * xe, 1.2 * xe], y=[0, 0], mode="lines",
                             name="Rotation axis", line=dict(color="black", width=1, dash="dashdot"),
                             hoverinfo="skip"))

    if dimensions:
        m = 0.10 * d2
        _extension(fig, (0, r2), (0, r2 + 1.15 * m))
        _extension(fig, (zE, r2), (zE, r2 + 1.15 * m))
        _extension(fig, (xe, r2), (xe, r2 + 0.65 * m))
        _dim(fig, (0, r2 + m), (zE, r2 + m), f"z<sub>E</sub> = {zE:g}")
        _dim(fig, (zE, r2 + 0.5 * m), (xe, r2 + 0.5 * m), f"b<sub>2</sub> = {b2:g}")

        _dim(fig, (0, r1 - 0.25 * m), (g1, r1 - 0.25 * m), "", text_pos=(g1 / 2, r1 - 0.55 * m))
        fig.add_annotation(x=g1 / 2, y=r1 - 0.55 * m, text=f"g<sub>1</sub> = {g1:g}",
                           showarrow=False, font=dict(color="#1f77b4", size=13))
        _extension(fig, (g1, r1), (g1, r1 - 0.35 * m))

        _extension(fig, (0, r1), (-1.2 * m, r1))
        _extension(fig, (xe, r2), (xe + 1.2 * m, r2))
        _dim(fig, (-m, 0), (-m, r1), f"d<sub>1</sub> = {d1:g}", text_pos=(-m, r1 / 2))
        _dim(fig, (xe + m, 0), (xe + m, r2), f"d<sub>2</sub> = {d2:g}", text_pos=(xe + m, r2 / 2))

        # shroud outlet angle: show the radial reference and the tangent
        L = 0.25 * m * 4
        e = np.radians(eps_ds)
        fig.add_trace(go.Scatter(
            x=[zE, zE, None, zE, zE + L * np.sin(e)], y=[r2, r2 + L, None, r2, r2 + L * np.cos(e)],
            mode="lines", line=dict(color="gray", width=1, dash="dot"), hoverinfo="skip",
            showlegend=False))
        fig.add_annotation(x=zE + 0.5 * L * np.sin(e) * 1.0, y=r2 + L * 1.05,
                           text=f"ε<sub>DS</sub> = {eps_ds:g}°", showarrow=False,
                           xanchor="left", font=dict(size=12, color="gray"))

    fig.update_layout(
        template=template,
        title=f"Impeller meridional cross-section (d2 = {d2:g}, b2 = {b2:g}, zE = {zE:g}, d1 = {d1:g})",
        xaxis=dict(title="axial direction", zeroline=False, constrain="domain"),
        yaxis=dict(title="radius  r = d/2", scaleanchor="x", scaleratio=1, zeroline=False),
        legend=dict(orientation="h", yanchor="bottom", y=1.02, x=0),
        hovermode="closest", height=800,
    )
    return fig


def plot_impeller_3d(d2, b2, zE, d1, g1, eps_ds, R_ds,
                     eps_ts=None, R_ts=None, sweep_deg=270, n_theta=90, template="plotly_white"):
    """Shroud and hub surfaces of revolution (cut-away by default)."""
    shroud, hub, _, _ = build_geometry(d2, b2, zE, d1, g1, eps_ds, R_ds, eps_ts, R_ts)
    th = np.radians(np.linspace(0, sweep_deg, n_theta))
    fig = go.Figure()
    for xy, name, col in ((shroud, "Shroud (DS)", "#d62728"), (hub, "Hub (TS)", "#1f77b4")):
        X = np.tile(xy[:, 0], (n_theta, 1))
        Y = np.outer(np.cos(th), xy[:, 1])
        Z = np.outer(np.sin(th), xy[:, 1])
        fig.add_trace(go.Surface(x=X, y=Y, z=Z, name=name, showscale=False, opacity=0.85,
                                 colorscale=[[0, col], [1, col]], showlegend=True))
    fig.update_layout(template=template, height=800, title="Surfaces of revolution",
                      scene=dict(aspectmode="data", xaxis_title="axial",
                                 yaxis_title="y", zaxis_title="z"))
    return fig


if __name__ == "__main__":
    params = dict(d2=400, b2=60, zE=110, d1=200, g1=20, eps_ds=5, R_ds=80)

    fig = plot_impeller(**params)
    #fig.write_html("impeller_section.html")      # shareable, interactive
    # fig.write_image("impeller_section.png", scale=2)   # needs kaleido
    fig.show()

    # plot_impeller_3d(**params).show()