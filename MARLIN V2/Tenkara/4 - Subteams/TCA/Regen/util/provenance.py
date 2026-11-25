"""
Records what made each output, so a file in output/ or a saved plot can be
traced back to its inputs and code.

One stamp per run (computed once, shared by every output of that run):
    time      local date and time of the run
    git       commit the code was run from, "+dirty" if there were uncommitted
              changes to tracked files (untracked files don't count)
    config    path of the yaml that was loaded and the first 12 hex digits of its
              SHA-256. Same hash = same inputs, byte for byte.
    script    the script that was run

Where it goes:
    output/inputs_used.yaml   exact copy of the loaded yaml, with the stamp and the
                              list of files this run wrote on top. Files in output/
                              that aren't on that list are left over from an
                              earlier run.
    contour.dxf               stamp in the DXF header (custom properties, not geometry)
    channel_heights.csv       stamp as '#' comment lines above the column header
    plots                     stamp as a small footer line
    .dat NX point files       nothing inline (NX reads every line as a point);
                              covered by inputs_used.yaml
"""
import datetime
import pathlib
import subprocess
import sys

from util import config as c

_written = []      # files written by this run, for inputs_used.yaml
_stamp = None


def _git_version():
    try:
        out = subprocess.run(["git", "-C", str(c.ROOT), "describe", "--always", "--dirty=+dirty"],
                             capture_output=True, text=True, timeout=10)
        return out.stdout.strip() if out.returncode == 0 else "no git"
    except (OSError, subprocess.SubprocessError):
        return "no git"


def stamp():
    """The run's stamp as a dict (time, git, config, config_sha256, script)."""
    global _stamp
    if _stamp is None:
        try:
            cfg_path = c.PATH.resolve().relative_to(c.ROOT).as_posix()
        except ValueError:
            cfg_path = str(c.PATH.resolve())
        try:
            script = pathlib.Path(sys.argv[0]).resolve().relative_to(c.ROOT).as_posix()
        except (ValueError, OSError):
            script = sys.argv[0]
        _stamp = dict(
            time=datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            git=_git_version(),
            config=cfg_path,
            config_sha256=c.SHA256[:12],
            script=script,
        )
    return _stamp


def stamp_line():
    """One line, for plot footers and file headers."""
    s = stamp()
    return (f"regen {s['time']} | git {s['git']} | {s['config']} sha256 {s['config_sha256']}"
            f" | {s['script']}")


def record(path):
    """Note that this run wrote `path` (listed in inputs_used.yaml)."""
    _written.append(pathlib.Path(path))


def stamp_figure(fig):
    """Put the stamp in a small grey footer, so it's in any image saved from the window."""
    fig.text(0.005, 0.003, stamp_line(), fontsize=6, color="0.45",
             ha="left", va="bottom")


def write_inputs_used(out_dir):
    """Write out_dir/inputs_used.yaml: stamp + files written + the exact yaml loaded."""
    out_dir = pathlib.Path(out_dir)
    out_dir.mkdir(exist_ok=True)
    s = stamp()
    lines = ["# Inputs used for the files listed below. Written by util/provenance.py.",
             f"#   time     {s['time']}",
             f"#   git      {s['git']}",
             f"#   config   {s['config']}  (sha256 {c.SHA256})",
             f"#   script   {s['script']}",
             "#   files written by this run:"]
    for p in _written:
        try:
            p = p.resolve().relative_to(c.ROOT).as_posix()
        except ValueError:
            pass
        lines.append(f"#     {p}")
    lines += ["# Anything else in this folder is from an earlier run.",
              "# ---------------- exact copy of the loaded yaml below ----------------"]
    path = out_dir / "inputs_used.yaml"
    with open(path, "w", encoding="utf-8", newline="") as f:     # newline="": copy stays byte-exact
        f.write("\n".join(lines) + "\n" + c.TEXT)
    print(f"  Inputs used written to {path}")
