"""Desktop GUI for CEA-based rocket chamber gas-side thermal analysis.

Install: pip install numpy pandas scipy matplotlib pint pyyaml
Also install the `cea` package used by your existing solver.
Run: python thermal_gui.py

The GUI reads an existing Inputs/TCA_params.yaml (or another selected YAML),
accepts x_mm/y_mm geometry or a precomputed properties.csv, plus a GUI wall-temperature input, and produces:
Outputs/properties.csv (CEA mode only), heat_transfer_parameters.csv,
ansys_input2.csv (headerless, 500 points). The contour input is never overwritten.
"""
from __future__ import annotations

import threading
import queue
import traceback
import os
from datetime import datetime
from pathlib import Path
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import numpy as np
import pandas as pd
import yaml
from scipy.interpolate import CubicSpline
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

COLUMNS = ('x_mm', 'y_mm')
PROJECT_ROOT = Path(__file__).resolve().parent.parent



def validate_geometry(df: pd.DataFrame) -> pd.DataFrame:
    if not set(COLUMNS).issubset(df.columns):
        raise ValueError('Geometry requires columns: x_mm, y_mm.')
    data = df[list(COLUMNS)].copy()
    for name in COLUMNS:
        data[name] = pd.to_numeric(data[name], errors='raise')
    data = data.sort_values('x_mm').reset_index(drop=True)
    if len(data) < 4:
        raise ValueError('Provide at least four geometry stations for a cubic spline.')
    if not np.isfinite(data.to_numpy()).all():
        raise ValueError('Geometry values must be finite.')
    if data.x_mm.duplicated().any():
        raise ValueError('x_mm must be unique at every station.')
    if (data.y_mm <= 0).any():
        raise ValueError('Radius y_mm must be positive.')
    if not (data.x_mm.min() <= 0 <= data.x_mm.max()):
        raise ValueError('The x range must include x=0 (the assumed throat location).')
    return data


def run_analysis(geometry, params_path, output_dir, wall_temperature, log, existing_properties=None,
                 chamber_roughness_m=6.35e-6):
    """Orchestrate the user's CEA and Bartz modules, without copying equations."""
    from properties_krishupdate import generate_properties, validate_geometry as check_geometry
    from bartz_krishupdate import calculate_bartz

    if existing_properties is None:
        g = check_geometry(geometry)
        log('Running properties_krishupdate.generate_properties ...')
        props = generate_properties(g, params_path, log=log)
        log('Running bartz_krishupdate.calculate_bartz ...')
        heat, props, ansys = calculate_bartz(
            props, wall_temperature, params_path, output_dir, log, contour_geometry=g,
            chamber_roughness_m=chamber_roughness_m)
    else:
        log('Using existing properties CSV; CEA skipped ...')
        # Do not overwrite the user's source properties file.
        heat, props, ansys = calculate_bartz(
            existing_properties, wall_temperature, params_path, output_dir, log,
            chamber_roughness_m=chamber_roughness_m)
        g = pd.DataFrame({'x_mm': props['x_m'] * 1000,
                          'y_mm': props['y_m'] * 1000})
        # Existing-properties mode reads the source CSV without overwriting it.
        # Geometry is used for preview only; do not export or modify contour.csv.
    return g, heat, props, ansys


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title('Rocket Chamber | CEA & Gas-Side Thermal Analysis')
        self.geometry('1230x770')
        self.minsize(950, 600)
        self.events = queue.Queue()
        self.last_result = None
        self.last_mode = None
        self.geometry_df = pd.DataFrame(columns=COLUMNS)
        self.mode = tk.StringVar(value='properties')
        self.properties_path = tk.StringVar(value='Outputs/properties.csv')
        self.wall_temp = tk.StringVar(value='800')
        self.chamber_roughness = tk.StringVar(value='6.35')
        self.lookup_x = tk.StringVar()
        self.lookup_result = tk.StringVar(value='Run an analysis to enable axial-position lookup.')
        self.yaml_path = tk.StringVar(value='Inputs/TCA_params.yaml')
        self.output_path = tk.StringVar(value='Outputs')
        self.status = tk.StringVar(value='Choose an existing properties CSV or import contour geometry.')
        # Make relative Inputs/Outputs paths resolve from the TCA root.
        os.chdir(PROJECT_ROOT)
        self._build_ui()
        self.after(100, self._poll)

    def _build_ui(self):
        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill='both', expand=True, padx=12, pady=8)
        self.setpoint_view = ttk.Frame(self.notebook, padding=12)
        root = ttk.Frame(self.notebook, padding=12)
        cooling = ttk.Frame(self.notebook, padding=12)
        self.bartz_view = root
        self.coolant_view = cooling
        self.notebook.add(self.setpoint_view, text='Setpoint Chooser')
        self.notebook.add(root, text='Bartz')
        self.notebook.add(cooling, text='Coolant / IPA')
        self._build_setpoint_ui(self.setpoint_view)
        self._build_cooling_ui(cooling)
        ttk.Label(root, text='Thrust Chamber Thermal Analysis', font=('Segoe UI', 18, 'bold')).pack(anchor='w', pady=(0, 10))
        paths = ttk.LabelFrame(root, text='Inputs and outputs', padding=8)
        paths.pack(fill='x')
        paths.columnconfigure(1, weight=1)
        ttk.Label(paths, text='Engine YAML').grid(row=0, column=0, sticky='w')
        ttk.Entry(paths, textvariable=self.yaml_path).grid(row=0, column=1, sticky='ew', padx=8)
        ttk.Button(paths, text='Browse', command=self._browse_yaml).grid(row=0, column=2)
        self.properties_widgets = (
            ttk.Label(paths, text='Properties CSV'),
            ttk.Entry(paths, textvariable=self.properties_path),
            ttk.Button(paths, text='Browse', command=self._browse_properties),
        )
        self.properties_widgets[0].grid(row=1, column=0, sticky='w')
        self.properties_widgets[1].grid(row=1, column=1, sticky='ew', padx=8)
        self.properties_widgets[2].grid(row=1, column=2)
        ttk.Label(paths, text='Wall temperature [K]').grid(row=2, column=0, sticky='w')
        ttk.Entry(paths, textvariable=self.wall_temp, width=16).grid(row=2, column=1, sticky='w', padx=8)
        ttk.Label(paths, text='Constant value applied to every station').grid(row=2, column=2, sticky='w')
        modes = ttk.Frame(paths)
        modes.grid(row=3, column=0, columnspan=3, sticky='w', pady=4)
        ttk.Radiobutton(modes, text='Use existing properties CSV (skip CEA)', variable=self.mode, value='properties', command=self._mode_changed).pack(side='left', padx=(0, 18))
        ttk.Radiobutton(modes, text='Generate properties from geometry (run CEA)', variable=self.mode, value='geometry', command=self._mode_changed).pack(side='left')
        ttk.Label(paths, text='Output folder').grid(row=4, column=0, sticky='w')
        ttk.Entry(paths, textvariable=self.output_path).grid(row=4, column=1, sticky='ew', padx=8)
        ttk.Button(paths, text='Browse', command=self._browse_output).grid(row=4, column=2)
        ttk.Label(paths, text='Chamber roughness [µm]').grid(row=5, column=0, sticky='w')
        ttk.Entry(paths, textvariable=self.chamber_roughness, width=16).grid(row=5, column=1, sticky='w', padx=8)
        ttk.Label(paths, text='Gnielinski upstream; default 250 microinch from Constant_BCs').grid(row=5, column=2, sticky='w')

        # Keep the primary action visible above the expanding geometry/plot region.
        action_bar = ttk.Frame(root)
        action_bar.pack(fill='x', pady=(10, 4))
        self.run_button = ttk.Button(action_bar, text='Run Thermal Analysis', command=self._run)
        self.run_button.pack(side='left')
        self.export_button = ttk.Button(action_bar, text='Export Results (timestamped)', command=self._export_results, state='disabled')
        self.export_button.pack(side='left', padx=(12, 0))
        self.clear_button = ttk.Button(action_bar, text='Clear', command=self._clear)
        self.clear_button.pack(side='left', padx=(12, 0))
        ttk.Label(action_bar, textvariable=self.status).pack(side='left', padx=16)


        # Query the last successfully calculated thermal data without saving files.
        lookup = ttk.LabelFrame(root, text='Axial-position lookup (Gnielinski chamber / Bartz nozzle)', padding=8)
        lookup.pack(fill='x', pady=(4, 2))
        ttk.Label(lookup, text='Axial position x [mm]').pack(side='left')
        entry = ttk.Entry(lookup, textvariable=self.lookup_x, width=15)
        entry.pack(side='left', padx=(8, 6))
        entry.bind('<Return>', lambda _event: self._lookup_position())
        self.lookup_button = ttk.Button(lookup, text='Look up', command=self._lookup_position, state='disabled')
        self.lookup_button.pack(side='left', padx=(0, 12))
        ttk.Label(lookup, textvariable=self.lookup_result).pack(side='left')

        panels = ttk.Panedwindow(root, orient='horizontal')
        panels.pack(fill='both', expand=True, pady=10)
        left = ttk.Frame(panels, padding=(0, 5))
        right = ttk.Frame(panels, padding=(5, 5))
        panels.add(left, weight=2)
        panels.add(right, weight=3)
        gbox = ttk.LabelFrame(left, text='Geometry stations (x and radius y in mm) — geometry/CEA mode', padding=8)
        gbox.pack(fill='both', expand=True)
        controls = ttk.Frame(gbox)
        controls.pack(fill='x', pady=(0, 8))
        for name, cmd in [('Import CSV', self._import), ('Add row', self._add), ('Edit row', self._edit), ('Delete', self._delete), ('Export CSV', self._export)]:
            ttk.Button(controls, text=name, command=cmd).pack(side='left', padx=(0, 4))
        frame = ttk.Frame(gbox)
        frame.pack(fill='both', expand=True)
        self.tree = ttk.Treeview(frame, columns=COLUMNS, show='headings', selectmode='extended')
        for col in COLUMNS:
            self.tree.heading(col, text=col)
            self.tree.column(col, width=115, anchor='center')
        scroll = ttk.Scrollbar(frame, command=self.tree.yview)
        self.tree.configure(yscrollcommand=scroll.set)
        self.tree.pack(side='left', fill='both', expand=True)
        scroll.pack(side='right', fill='y')
        self.tree.bind('<Double-1>', lambda _event: self._edit())
        ttk.Label(gbox, text='Geometry mode: points on both sides of x=0; x_mm/y_mm or x_m/y_m accepted.').pack(anchor='w', pady=(7, 0))
        self.figure = Figure(figsize=(7, 6), dpi=100)
        self.ax = self.figure.add_subplot(111)
        self.ax.set_title('Contour preview')
        self.ax.set_xlabel('Axial position x [mm]')
        self.ax.set_ylabel('Radius y [mm]')
        # One mm horizontally occupies the same screen distance as one mm vertically.
        self.ax.set_aspect('equal', adjustable='box')
        self.ax.grid(alpha=0.25)
        self.canvas = FigureCanvasTkAgg(self.figure, master=right)
        self.canvas.get_tk_widget().pack(fill='both', expand=True)
        logbox = ttk.LabelFrame(root, text='Run log', padding=4)
        logbox.pack(fill='x', pady=(8, 0))
        self.log = tk.Text(logbox, height=5, state='disabled', wrap='word')
        self.log.pack(fill='x')

    def _build_setpoint_ui(self, root):
        ttk.Label(root, text='Setpoint Chooser', font=('Segoe UI', 18, 'bold')).pack(anchor='w')
        ttk.Label(root, text='Setpoint selection controls will be added here.').pack(anchor='w', pady=6)

    def _build_cooling_ui(self, root):
        self.cooling_result = None
        self.cooling_inputs = {}
        ttk.Label(root, text='IPA Cooling-Channel Heat Transfer', font=('Segoe UI', 18, 'bold')).pack(anchor='w')
        ttk.Label(root, text='Constant coolant properties and geometry; equal flow through parallel channels.').pack(anchor='w', pady=6)
        panels = ttk.Panedwindow(root, orient='horizontal')
        panels.pack(fill='both', expand=True)
        left = ttk.Frame(panels)
        right = ttk.Frame(panels)
        panels.add(left, weight=2)
        panels.add(right, weight=3)
        groups = [
            ('Channel geometry (blank until supplied)', [
                ('area', 'Flow area per channel [mm²]'),
                ('perimeter', 'Wetted perimeter per channel [mm]'),
                ('channels', 'Number of parallel channels'),
                ('roughness', 'Absolute roughness [mm]'),
                ('x_start', 'Axial start x [mm]'), ('x_end', 'Axial end x [mm]')]),
            ('IPA bulk properties at operating conditions (manual)', [
                ('mdot', 'Total coolant mass flow [kg/s]'),
                ('rho', 'Density [kg/m³]'), ('mu', 'Dynamic viscosity [Pa·s]'),
                ('cp', 'Specific heat [J/(kg·K)]'),
                ('k', 'Thermal conductivity [W/(m·K)]')])]
        for title, fields in groups:
            box = ttk.LabelFrame(left, text=title, padding=8)
            box.pack(fill='x', pady=6)
            for row, (key, label) in enumerate(fields):
                value = tk.StringVar()
                self.cooling_inputs[key] = value
                ttk.Label(box, text=label).grid(row=row, column=0, sticky='w', pady=3)
                ttk.Entry(box, textvariable=value, width=16).grid(row=row, column=1, padx=8)
        actions = ttk.Frame(root)
        actions.pack(fill='x', pady=8)
        ttk.Button(actions, text='Run Cooling Analysis', command=self._run_cooling).pack(side='left')
        self.cooling_export = ttk.Button(actions, text='Export Results (timestamped)', command=self._export_cooling, state='disabled')
        self.cooling_export.pack(side='left', padx=8)
        ttk.Button(actions, text='Clear', command=self._clear_cooling).pack(side='left')
        lookup = ttk.LabelFrame(root, text='Axial-position lookup', padding=8)
        lookup.pack(fill='x')
        self.cooling_x = tk.StringVar()
        ttk.Label(lookup, text='Axial position x [mm]').pack(side='left')
        ttk.Entry(lookup, textvariable=self.cooling_x, width=14).pack(side='left', padx=8)
        ttk.Button(lookup, text='Look up', command=self._lookup_cooling).pack(side='left')
        self.cooling_summary = tk.StringVar(value='Enter channel geometry and IPA properties to calculate HTC.')
        ttk.Label(root, textvariable=self.cooling_summary, wraplength=1100).pack(anchor='w', pady=8)
        self.cooling_figure = Figure(figsize=(7, 5), dpi=100)
        self.cooling_ax = self.cooling_figure.add_subplot(111)
        self.cooling_canvas = FigureCanvasTkAgg(self.cooling_figure, master=right)
        self.cooling_canvas.get_tk_widget().pack(fill='both', expand=True)
        self._plot_cooling()

    def _plot_cooling(self):
        ax = self.cooling_ax
        ax.clear()
        ax.set(title='IPA coolant convective heat transfer', xlabel='Axial position x [mm]',
               ylabel='Convective HTC h_c [W/(m²·K)]')
        if self.cooling_result is not None:
            for column, label in [('h_gnielinski_W_m2K', 'Gnielinski'), ('h_dittus_boelter_W_m2K', 'Dittus–Boelter')]:
                if self.cooling_result[column].notna().any():
                    ax.plot(self.cooling_result.x_mm, self.cooling_result[column], label=label,
                            linestyle='--' if 'dittus' in column else '-')
            ax.legend()
        ax.grid(alpha=.25)
        ax.ticklabel_format(axis='y', style='plain', useOffset=False)
        self.cooling_figure.tight_layout()
        self.cooling_canvas.draw_idle()

    def _run_cooling(self):
        from constant_bc_loader import load_correlations
        import contextlib
        import io
        self.cooling_result = None
        self.cooling_export.configure(state='disabled')
        self._plot_cooling()
        try:
            v = {key: float(value.get()) for key, value in self.cooling_inputs.items()}
            if not all(np.isfinite(value) for value in v.values()):
                raise ValueError('All inputs must be finite.')
            for key in ('area', 'perimeter', 'channels', 'mdot', 'rho', 'mu', 'cp', 'k'):
                if v[key] <= 0:
                    raise ValueError(f'{key} must be positive.')
            if not v['channels'].is_integer() or v['roughness'] < 0 or v['x_end'] <= v['x_start']:
                raise ValueError('Use an integer channel count, nonnegative roughness, and axial end greater than start.')
            bc = load_correlations()
            area = v['area'] * 1e-6
            diameter = 4 * area / (v['perimeter'] * 1e-3)
            velocity = v['mdot'] / (v['channels'] * v['rho'] * area)
            re = v['rho'] * velocity * diameter / v['mu']
            pr = v['mu'] * v['cp'] / v['k']
            warnings = io.StringIO()
            with contextlib.redirect_stdout(warnings):
                friction = bc.friction_factor(v['roughness'] * 1e-3, diameter, re)
                nu_g = bc.gnielinski(friction, re, pr) if friction is not None else None
                nu_d = bc.dittus_boelter(re, pr, fluid_being_cooled=False)
            messages = warnings.getvalue().strip()
            if nu_g is None and nu_d is None:
                raise ValueError('Neither correlation is valid for these inputs.\n' + messages)
            x = np.linspace(v['x_start'], v['x_end'], 500)
            self.cooling_result = pd.DataFrame({
                'x_mm': x, 'h_gnielinski_W_m2K': nu_g * v['k'] / diameter if nu_g is not None else np.nan,
                'h_dittus_boelter_W_m2K': nu_d * v['k'] / diameter if nu_d is not None else np.nan,
                'Re': re, 'Pr': pr, 'Darcy_f': friction if friction is not None else np.nan,
                'Nu_gnielinski': nu_g if nu_g is not None else np.nan,
                'Nu_dittus_boelter': nu_d if nu_d is not None else np.nan,
                'hydraulic_diameter_m': diameter, 'velocity_m_s': velocity, **v})
            self.cooling_summary.set(f'D_h = {diameter*1000:.4g} mm; velocity = {velocity:.4g} m/s; Re = {re:,.0f}; Pr = {pr:.4g}. '
                                     'Constant conditions give horizontal HTC curves. ' + messages)
            self.cooling_export.configure(state='normal')
            self._plot_cooling()
        except Exception as exc:
            self.cooling_summary.set('Calculation failed; check inputs and Constant_BCs dependencies.')
            messagebox.showerror('Cooling analysis failed', str(exc))

    def _lookup_cooling(self):
        try:
            if self.cooling_result is None:
                raise ValueError('Run a cooling analysis first.')
            x = float(self.cooling_x.get())
            data = self.cooling_result
            if not np.isfinite(x) or not data.x_mm.iloc[0] <= x <= data.x_mm.iloc[-1]:
                raise ValueError('Position must lie within the calculated axial range.')
            values = []
            for column, label in [('h_gnielinski_W_m2K', 'Gnielinski'), ('h_dittus_boelter_W_m2K', 'Dittus–Boelter')]:
                h = float(data[column].iloc[0])
                values.append(f'{label}: {h:,.2f} W/(m²·K)' if np.isfinite(h) else f'{label}: outside validity range')
            self.cooling_summary.set(f'x = {x:g} mm; ' + '; '.join(values))
        except ValueError as exc:
            self.cooling_summary.set(str(exc))

    def _export_cooling(self):
        if self.cooling_result is None:
            return
        folder = filedialog.askdirectory(title='Export coolant HTC results')
        if not folder:
            return
        try:
            target = Path(folder) / f'cooling_htc_{datetime.now():%Y%m%d_%H%M%S_%f}.csv'
            with target.open('x', newline='', encoding='utf-8') as fh:
                self.cooling_result.to_csv(fh, index=False)
            self.cooling_summary.set(f'Exported {target}')
        except Exception as exc:
            messagebox.showerror('Export failed', str(exc))

    def _clear_cooling(self):
        for value in self.cooling_inputs.values():
            value.set('')
        self.cooling_x.set('')
        self.cooling_result = None
        self.cooling_export.configure(state='disabled')
        self.cooling_summary.set('Enter channel geometry and IPA properties to calculate HTC.')
        self._plot_cooling()

    def _clear(self):
        """Restore the startup inputs, empty geometry, and unused analysis state."""
        self.mode.set('properties')
        self._update_properties_visibility()
        self.properties_path.set('Outputs/properties.csv')
        self.wall_temp.set('800')
        self.chamber_roughness.set('6.35')
        self.yaml_path.set('Inputs/TCA_params.yaml')
        self.output_path.set('Outputs')
        self.geometry_df = pd.DataFrame(columns=COLUMNS)
        self.last_result = None
        self.last_mode = None
        self.export_button.configure(state='disabled')
        self.lookup_button.configure(state='disabled')
        self.lookup_x.set('')
        self.lookup_result.set('Run an analysis to enable axial-position lookup.')
        self.log.configure(state='normal')
        self.log.delete('1.0', 'end')
        self.log.configure(state='disabled')
        self.run_button.configure(state='normal')
        self._refresh()
        self.status.set('Choose an existing properties CSV or import contour geometry.')

    def _browse_properties(self):
        val = filedialog.askopenfilename(filetypes=[('CSV', '*.csv')])
        if val:
            self.properties_path.set(val)

    def _update_properties_visibility(self):
        for widget in self.properties_widgets:
            if self.mode.get() == 'properties':
                widget.grid()
            else:
                widget.grid_remove()

    def _mode_changed(self):
        self._update_properties_visibility()
        self.last_result = None
        self.export_button.configure(state='disabled')
        self.lookup_button.configure(state='disabled')
        self.lookup_result.set('Run an analysis to enable axial-position lookup.')
        self.status.set('CEA will be skipped.' if self.mode.get() == 'properties' else 'Geometry mode: import x_mm/y_mm or x_m/y_m CSV, or add rows.')

    def _browse_yaml(self):
        val = filedialog.askopenfilename(filetypes=[('YAML', '*.yaml *.yml'), ('All files', '*')])
        if val:
            self.yaml_path.set(val)

    def _browse_output(self):
        val = filedialog.askdirectory()
        if val:
            self.output_path.set(val)

    def _import(self):
        val = filedialog.askopenfilename(filetypes=[('CSV', '*.csv')])
        if not val:
            return
        try:
            from properties_krishupdate import read_geometry_csv
            self.geometry_df = read_geometry_csv(val)
            self._refresh()
            self.status.set(f'Imported {len(self.geometry_df)} stations')
        except Exception as exc:
            messagebox.showerror('CSV import failed', str(exc))

    def _row_dialog(self, title, defaults=None):
        dialog = tk.Toplevel(self)
        dialog.title(title)
        dialog.resizable(False, False)
        dialog.transient(self)
        dialog.grab_set()
        values = []
        for i, col in enumerate(COLUMNS):
            ttk.Label(dialog, text=col).grid(row=i, column=0, padx=12, pady=7)
            field = ttk.Entry(dialog, width=25)
            field.grid(row=i, column=1, padx=12, pady=7)
            if defaults is not None:
                field.insert(0, str(defaults[i]))
            values.append(field)
        result = []
        def save():
            try:
                result.extend(float(v.get()) for v in values)
                if not np.isfinite(result).all():
                    raise ValueError('All values must be finite.')
                dialog.destroy()
            except ValueError:
                result.clear()
                messagebox.showerror('Invalid row', 'Enter finite numeric values.', parent=dialog)
        ttk.Button(dialog, text='Save', command=save).grid(row=2, column=0, columnspan=2, pady=10)
        values[0].focus_set()
        self.wait_window(dialog)
        return result if len(result) == 2 else None

    def _add(self):
        row = self._row_dialog('Add geometry station')
        if row is not None:
            self.geometry_df.loc[len(self.geometry_df)] = row
            self._refresh()

    def _edit(self):
        selection = self.tree.selection()
        if not selection:
            return
        idx = int(selection[0])
        row = self._row_dialog('Edit geometry station', self.geometry_df.iloc[idx].tolist())
        if row is not None:
            self.geometry_df.iloc[idx] = row
            self._refresh()

    def _delete(self):
        ids = [int(item) for item in self.tree.selection()]
        if ids:
            self.geometry_df = self.geometry_df.drop(self.geometry_df.index[ids]).reset_index(drop=True)
            self._refresh()

    def _export(self):
        """Explicitly export a timestamped geometry copy, without overwrite."""
        if self.geometry_df.empty:
            messagebox.showwarning('No geometry', 'Import or enter geometry first.')
            return
        folder = filedialog.askdirectory(title='Choose folder for geometry export')
        if not folder:
            return
        stamp = datetime.now().strftime('%Y%m%d_%H%M%S_%f')
        target = Path(folder) / f'geometry_{stamp}.csv'
        try:
            with target.open('x', newline='', encoding='utf-8') as fh:
                self.geometry_df.to_csv(fh, index=False)
            messagebox.showinfo('Geometry exported', str(target))
        except Exception as exc:
            messagebox.showerror('Export failed', str(exc))

    def _refresh(self):
        self.tree.delete(*self.tree.get_children())
        for i, row in self.geometry_df.reset_index(drop=True).iterrows():
            self.tree.insert('', 'end', iid=str(i), values=[f'{float(v):.6g}' for v in row])
        # Recreate the axes so no thermal twin-axis state survives a reset
        # or a new contour import.
        self.figure.clear()
        self.ax = self.figure.add_subplot(111)
        if len(self.geometry_df):
            try:
                g = self.geometry_df.astype(float).sort_values('x_mm')
                self.ax.plot(g.x_mm, g.y_mm, 'o-', label='Radius')
            except ValueError:
                pass
        self.ax.set(title='Contour preview', xlabel='Axial position x [mm]', ylabel='Radius y [mm]')
        # True 1:1 geometry scale (equal x and y lengths on the screen).
        # This setting applies only to contour previews, not h_g/T_r plots.
        self.ax.set_aspect('equal', adjustable='box')
        self.ax.grid(alpha=.25)
        self.figure.tight_layout()
        self.canvas.draw_idle()

    def _lookup_position(self):
        """Interpolate h_g, T_r and radius at a user-selected x in mm."""
        if self.last_result is None:
            self.lookup_result.set('Run an analysis first.')
            return
        try:
            x_mm = float(self.lookup_x.get())
            if not np.isfinite(x_mm):
                raise ValueError('Position must be finite.')
            x_query = x_mm / 1000.0
            _, heat, props, _ = self.last_result
            # Bartz returns props restricted to the active nozzle domain.
            h = heat.sort_values('x_m')
            g = props.sort_values('x_m')
            xx = h['x_m'].to_numpy(dtype=float)
            gx = g['x_m'].to_numpy(dtype=float)
            if x_query < xx[0] or x_query > xx[-1]:
                self.lookup_result.set(
                    f'Outside calculated region: {xx[0]*1000:.3f} to {xx[-1]*1000:.3f} mm.'
                )
                return
            # Piecewise-linear interpolation avoids cubic-spline overshoot.
            transition = h.loc[h['correlation'] == 'Bartz', 'x_m'].min()
            model = 'Gnielinski' if x_query < transition else 'Bartz'
            region = h.loc[h['correlation'] == model]
            hg = np.interp(x_query, region.x_m, region.hg_W_m2K)
            tr = np.interp(x_query, xx, h['Tr_K'].to_numpy(dtype=float))
            y_mm = np.interp(x_query, gx, g['y_m'].to_numpy(dtype=float)) * 1000.0
            self.lookup_result.set(
                f'{model} h_g: {hg:,.2f} W/m²·K     T_r: {tr:,.2f} K     y: {y_mm:.4f} mm'
            )
        except (ValueError, TypeError) as exc:
            self.lookup_result.set(f'Invalid axial position: {exc}')

    def _append_log(self, value):
        self.log.configure(state='normal')
        self.log.insert('end', value + '\n')
        self.log.see('end')
        self.log.configure(state='disabled')

    def _run(self):
        try:
            data = validate_geometry(self.geometry_df) if self.mode.get() == 'geometry' else None
            wall_temperature = float(self.wall_temp.get())
            roughness = float(self.chamber_roughness.get()) * 1e-6
            if not np.isfinite(roughness) or roughness < 0:
                raise ValueError('Chamber roughness must be finite and nonnegative.')
            if not np.isfinite(wall_temperature) or wall_temperature <= 0:
                raise ValueError('Enter a positive wall temperature in kelvin.')
            existing_properties = None
            if self.mode.get() == 'properties':
                existing_properties = Path(self.properties_path.get())
                if not existing_properties.is_file():
                    raise FileNotFoundError(f'Properties CSV not found: {existing_properties}')
            yaml_path = Path(self.yaml_path.get())
            if not yaml_path.is_file():
                raise FileNotFoundError(f'YAML file not found: {yaml_path}')
        except Exception as exc:
            messagebox.showerror('Invalid inputs', str(exc))
            return
        self.run_button.configure(state='disabled')
        self.clear_button.configure(state='disabled')
        self.export_button.configure(state='disabled')
        self.last_result = None
        self.lookup_button.configure(state='disabled')
        self.lookup_result.set('Calculating...')
        self.last_mode = self.mode.get()
        self.status.set('Running thermal analysis…')
        self._append_log('Starting analysis…')
        def work():
            try:
                result = run_analysis(data, str(yaml_path), self.output_path.get(), wall_temperature,
                                      lambda msg: self.events.put(('log', msg)),
                                      existing_properties=str(existing_properties) if existing_properties else None,
                                      chamber_roughness_m=roughness)
                self.events.put(('done', result))
            except Exception:
                self.events.put(('error', traceback.format_exc()))
        threading.Thread(target=work, daemon=True).start()

    def _export_results(self):
        """Save timestamped result files, never overwriting any existing path."""
        if self.last_result is None:
            messagebox.showwarning('No results', 'Run the calculation first.')
            return
        try:
            output = Path(self.output_path.get()).expanduser().resolve()
            # Saving only occurs after an explicit click on Export.
            output.mkdir(parents=True, exist_ok=True)
            _, heat, props, (x_ansys, hg_ansys, tr_ansys) = self.last_result
            timestamp = datetime.now().strftime('%Y%m%d_%H%M%S_%f')
            outputs = {
                f'heat_transfer_parameters_{timestamp}.csv': heat,
                f'ansys_input2_{timestamp}.csv': np.column_stack((x_ansys, hg_ansys, tr_ansys)),
            }
            if self.last_mode == 'geometry':
                outputs[f'properties_{timestamp}.csv'] = props
            # Exclusive creation guarantees that even coincident timestamps
            # cannot overwrite existing files.
            created = []
            for filename, data in outputs.items():
                target = output / filename
                with target.open('x', newline='', encoding='utf-8') as fh:
                    created.append(target)
                    if isinstance(data, pd.DataFrame):
                        data.to_csv(fh, index=False)
                    else:
                        np.savetxt(fh, data, delimiter=',')
            self.status.set(f'Exported {len(created)} timestamped files to {output}')
            self._append_log('Exported: ' + ', '.join(path.name for path in created))
            messagebox.showinfo('Export complete', '\n'.join(path.name for path in created))
        except Exception as exc:
            messagebox.showerror('Export failed', str(exc))

    def _poll(self):
        try:
            while True:
                kind, payload = self.events.get_nowait()
                if kind == 'log':
                    self._append_log(payload)
                elif kind == 'error':
                    self._append_log(payload)
                    self.status.set('Analysis failed')
                    self.lookup_result.set('Calculation failed; run again to look up values.')
                    self.run_button.configure(state='normal')
                    self.clear_button.configure(state='normal')
                    messagebox.showerror('Analysis failed', payload.splitlines()[-1])
                elif kind == 'done':
                    self.last_result = payload
                    _, heat, _, _ = payload
                    for extra in list(self.figure.axes):
                        if extra is not self.ax:
                            extra.remove()
                    self.ax.clear()
                    # The contour view uses equal x/y units. Thermal axes have
                    # unrelated units, so reset that geometry-only aspect rule.
                    # Otherwise Matplotlib zooms to an extremely narrow h_g band
                    # and the curve looks like two vertical spikes.
                    self.ax.set_aspect('auto')
                    ax2 = self.ax.twinx()
                    bartz = heat.loc[heat.correlation == 'Bartz']
                    chamber = heat.loc[heat.correlation == 'Gnielinski']
                    line_hg, = self.ax.plot(bartz.x_m, bartz.hg_W_m2K, label='Bartz h_g', color='tab:blue')
                    line_tr, = ax2.plot(heat.x_m, heat.Tr_K, label='Recovery temperature T_r [K]', color='tab:orange')
                    lines = [line_hg, line_tr]
                    if len(chamber):
                        line_chamber, = self.ax.plot(chamber.x_m, chamber.hg_W_m2K, label='Chamber Gnielinski h_g', color='tab:purple')
                        lines.insert(0, line_chamber)
                    if 'hg_custom_W_m2K' in heat:
                        line_custom, = self.ax.plot(
                            heat.x_m, heat.hg_custom_W_m2K,
                            label='RPE comparison h_g', color='tab:green', linestyle='--')
                        lines.insert(1, line_custom)
                    self.ax.legend(handles=lines, loc='best')
                    self.ax.set(xlabel='Axial position x [m]', ylabel='h_g [W/m²/K]', title='Thermal results')
                    ax2.set_ylabel('Recovery temperature [K]')
                    self.ax.grid(alpha=.25)
                    # Show actual h_g values instead of a misleading +1.6646e4
                    # offset with tiny decimal tick labels.
                    self.ax.ticklabel_format(axis='y', style='plain', useOffset=False)
                    self.figure.tight_layout()
                    self.canvas.draw_idle()
                    self.export_button.configure(state='normal')
                    self.lookup_button.configure(state='normal')
                    self.lookup_result.set('Enter an x position in mm and click Look up.')
                    self.status.set('Calculation complete — results plotted; no files saved')
                    self.run_button.configure(state='normal')
                    self.clear_button.configure(state='normal')
        except queue.Empty:
            pass
        self.after(100, self._poll)


if __name__ == '__main__':
    App().mainloop()
