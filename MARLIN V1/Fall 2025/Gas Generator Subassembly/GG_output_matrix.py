# Maelstrom Torch Testing Grapher
# Authors: Dominik Sloup
# First Created: 04/22/2023
# Last Updated: 04/25/2023
# Calculations done in SI units

import math
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import TextBox
import matlab.engine as m

try:
    m = m.connect_matlab()
except m.EngineError:
    print("Matlab engine not running. Starting a new session...")
    m = m.start_matlab()

m.addpath(r'C:\Users\miles\Onedrive\Documents\PURPL\cea_')
    
# ──────────────────────────────────────────────────────────────
#  CONSTANTS
# ──────────────────────────────────────────────────────────────
gamma_ox = 1.4
gamma_rp = 1.22
R = 8.314                               # J /(mol·K)
mm_ox, mm_rp = 32.00, 164.6            # g/mol
R_ox = R / (mm_ox/1000)                # J/(kg·K)
R_rp = R / (mm_rp/1000)              # J/(kg·K)
C_star     = 1056.8                   # m/s  (CEA)
T0         = 298.15                    # K
crit_ratio_ox = (2/(gamma_ox+1))**(gamma_ox/(gamma_ox-1))

relaxation = 0.05

psi_to_pa  = 6894.76

max_stiff = 0.2

cont_pc = None
cont_of = None
hat = None
pc_labels = []
of_labels = []
num_orifices = 8

# feed-pressure grid
Pox_psi = np.linspace(1, 1000, num=200)
Prp_psi = np.linspace(1, 1000, num=200)
P_ox_pa, P_rp_pa = Pox_psi*psi_to_pa, Prp_psi*psi_to_pa
X, Y = np.meshgrid(Pox_psi, Prp_psi)
rho_rp = 800 #kg/m^3

# choking coefficients
#K_rp = math.sqrt(gamma_rp/(R_rp*T0)) * (2/(gamma_rp+1))**((gamma_rp+1)/(2*(gamma_rp-1)))
Cd_rp = 0.7   # estimate
Cd_ox = 0.6   #estimate
Beta_rp = 0.07398/0.25    #d/D
#K_rp = Cd_rp/(np.sqrt(1-Beta_rp**4))
K_o = math.sqrt(gamma_ox/(R_ox*T0)) * (2/(gamma_ox+1))**((gamma_ox+1)/(2*(gamma_ox-1)))

# ──────────────────────────────────────────────────────────────
#  CORE CALCULATIONS
# ──────────────────────────────────────────────────────────────
def Mdot_rp(Prp, Pc, Afu):
    """Mass flow rate through RP-1 orifice."""
    mdot = np.zeros_like(Pc)
    mdot = np.where(
        (Prp > Pc) & (Prp - Pc != 0),
        Cd_rp * Afu * np.sqrt(2 * rho_rp * (Prp - Pc)),
        0
    )
    return np.maximum(mdot, 0)

def Mdot_ox(Pox, Pc, Ao):
    """Mass flow rate through oxidizer orifice."""
    mdot = np.zeros_like(Pc)
    rho_ox = (Pox / (R_ox * T0))
    mdot = np.where((Pc / Pox < crit_ratio_ox), 
            Cd_ox * Ao * np.sqrt(2 * rho_ox *(Pox) * (gamma_ox/(gamma_ox-1)) * ((Pc/Pox)**(2/gamma_ox) - (Pc/ Pox)**((gamma_ox+1)/gamma_ox))),
            Cd_ox * Ao * Pox * K_o
        )
        
    return np.maximum(mdot, 0)

def Find_Temp(Pc, OF):
    # Ensure arrays are float and contiguous
    print("Pc:", Pc.dtype)
    Pc1 = np.asarray(Pc, dtype=np.float64).tolist()
    OF1 = np.asarray(OF, dtype=np.float64).tolist()
    
    mPc = Pc1
    mOF = OF1
    # out = m.Temp_arr(mPc, mOF)
    # if out:
    #     print('worked')

def compute_fields(d_fuel, d_ox, d_exit):
    inch = 0.0254
    Afu = num_orifices * (math.pi * (d_fuel * inch / 2) ** 2)
    Ao = num_orifices * (math.pi * (d_ox * inch / 2) ** 2)
    Ae = math.pi * (d_exit * inch / 2) ** 2
    error = []

    # Prepare arrays
    Pc_pa = np.full((len(Prp_psi), len(Pox_psi)), 1 * psi_to_pa)
    mdot_ox = Cd_ox * Ao * P_ox_pa[None, :] * K_o

    for i in range(200):  # Iterate for convergence
        mdot_rp = Mdot_rp(P_rp_pa[:, None], Pc_pa, Afu)
        mdot_ox = Mdot_ox(P_ox_pa[None, :], Pc_pa, Ao)
                                  
        Pc_pa_new = (mdot_ox + mdot_rp) * (C_star / Ae)
        # Under-relaxation to help convergence
        Pc_pa = (1-relaxation) * Pc_pa + relaxation * Pc_pa_new
        error.append((np.abs(Pc_pa[100,160] - Pc_pa_new[100,160])))
    mask_error = np.abs(Pc_pa - Pc_pa_new) > 0.5
    # test = np.array(m.test(4))
    # print("test:", test)
    # Ensure Pc has the correct shape
    Pc = Pc_pa / psi_to_pa
    OF = np.where(
        Mdot_rp(P_rp_pa[:,None], Pc_pa, Afu) != 0,
        mdot_ox / Mdot_rp(P_rp_pa[:, None], Pc_pa, Afu), 
        100
        )
    vapor_p_rp = 2000 / psi_to_pa #psi
    rp_p_loss = (0.5*(rho_rp) * (mdot_rp/(rho_rp*Afu))**2)  / psi_to_pa   #rp-1 dynamic rpessure loss
    # plt.figure(1)
    # plt.plot(error)
    # plt.title('Convergence of Chamber Pressure')
    mask = (((X - Pc) / Pc) <= max_stiff) | (Y - rp_p_loss <= vapor_p_rp)  # unchoke threshold
    #mask = (Y - rp_p_loss <= vapor_p_rp)  # unchoke threshold
    # plt.figure(2)
    # plt.plot(mdot_ox)
    # plt.title('Mass Flow Rate of Oxidizer')
    return Pc, OF, mask, mask_error
    

def pick_levels(arr):
    return np.linspace(np.nanmin(arr), np.nanmax(arr), 5)[1:4]      # 25–75 %

# ──────────────────────────────────────────────────────────────
#  INITIAL VALUES
# ──────────────────────────────────────────────────────────────
d_fuel0, d_ox0, d_exit0 = 0.057245, 0.07398, 0.56536       # in         
Pc0, OF0, mask0, mask_error = compute_fields(d_fuel0, d_ox0, d_exit0)
# Find_Temp(Pc0, OF0)


# ──────────────────────────────────────────────────────────────
#  FIGURE & MAIN AXES  (≈1000 px wide, *not* maximised)
# ──────────────────────────────────────────────────────────────
fig = plt.figure(figsize=(10, 10), dpi=70)               # 10 in × 100 dpi ≈ 1000 px
main_ax = fig.add_axes([0.05, 0.42, 0.88, 0.50])          # [L, B, W, H]

pcm = main_ax.pcolormesh(X, Y, Pc0, shading="auto")
cbar = fig.colorbar(pcm, ax=main_ax, label="P₍c,eff₎ (psi)")

cont_pc   = main_ax.contour(X, Y, Pc0, levels = pick_levels(Pc0),
                            colors='k', linewidths=1)
pc_labels = main_ax.clabel(cont_pc, fmt='%d psi', inline=True, fontsize=14)

cont_of   = main_ax.contour(X, Y, OF0, levels = [0.4, 0.6, 0.8, 1.0],
                            colors='k', linewidths=1)
of_labels = main_ax.clabel(cont_of, fmt='OF: %.2f', inline=True, fontsize=14)

hat1 = main_ax.contourf(X, Y, mask0, levels=[0.5,1.0], colors=['none'],
                       hatches=['\\\\'], zorder=2, edgecolor='red')

hat2 = main_ax.contourf(X, Y, mask_error, levels=[0.5,1.0], colors=['none'], 
                          hatches=['\\\\'], zorder=3)
                    
#print("Unique values in mask0:", np.unique(mask0))
#print(hat.collections)
if hasattr(hat, 'collections'):
    for coll in hat.collections:
       coll.set(edgecolor='red', linewidth=0.5)

main_ax.set(xlim=(0,1000), ylim=(0,1000), aspect='equal',
            xlabel='Oxidizer Feed Pressure (psi)',
            ylabel='Fuel Feed Pressure (psi)',
            title='Injector Choking Map  &  Chamber-Pressure Contours')

point_handle = None     # placeholder for marker

# ──────────────────────────────────────────────────────────────
#  TEXT-BOX LAYOUT
# ──────────────────────────────────────────────────────────────
BOX_W, BOX_H = 0.18, 0.06
H_SP         = 0.07
ROW1_Y, ROW2_Y = 0.25, 0.10            # ← widened vertical gap


# horizontal centring under plot+colour-bar
fig.canvas.draw()
left  = min(main_ax.get_position().x0, cbar.ax.get_position().x0)
right = max(main_ax.get_position().x1, cbar.ax.get_position().x1)
X_START = (left + right - (3*BOX_W + 2*H_SP)) / 2

def make_box(x, y, label, init=''):
    ax_box = fig.add_axes([x, y, BOX_W, BOX_H])
    ax_box.text(0.5, 1.2, label, transform=ax_box.transAxes,
                ha='center', va='bottom')
    return TextBox(ax_box, '', initial=init)

# row 1 – diameters
tb_h = make_box(X_START,                  ROW1_Y, 'Fuel Orifice diameter (in)',  str(d_fuel0))
tb_o = make_box(X_START+BOX_W+H_SP,       ROW1_Y, 'Ox Orifice Diameter (in)',    str(d_ox0))
tb_e = make_box(X_START+2*(BOX_W+H_SP),   ROW1_Y, 'Throat Diameter (in)',  str(d_exit0))

# row 2 – pressures
pf_box = make_box(X_START,                  ROW2_Y, 'Fuel Line Pressure (psi)', '')
po_box = make_box(X_START+BOX_W+H_SP,       ROW2_Y, 'Oxidizer Line Pressure (psi)', '')
pc_box = make_box(X_START+2*(BOX_W+H_SP),   ROW2_Y, 'Chamber Pressure (psi)',  '')

# OF read-out (bottom-most)
of_text = fig.text(X_START+BOX_W,0.04,'OF = ', va='center', fontsize=12, weight='bold')

# ──────────────────────────────────────────────────────────────
#  CALLBACKS
# ──────────────────────────────────────────────────────────────
def compute_inlet(_):
    """Solve third pressure, update marker/OF, and display unchoke threshold."""
    global point_handle
    dh, dox, de = map(float, (tb_h.text, tb_o.text, tb_e.text))
    inch = 0.0254
    Ah, Ao, Ae = [math.pi*(d*inch/2)**2 for d in (dh, dox, de)]
    Ah = num_orifices * Ah
    Ao = num_orifices * Ao

    pf = float(pf_box.text) if pf_box.text else None
    po = float(po_box.text) if po_box.text else None
    pc = float(pc_box.text) if pc_box.text else None
    
    if pf is not None and po is not None:

        pc = Pc0[int(po/5), int(pf/5)]
        pc_box.set_val(f'{pc:.2f}')
        
    # elif pf is not None and pc is not None:
    #     po = (pc*(Ae/C_star) - Ah*K_rp*pf) / (Ao*K_o)
    #     po_box.set_val(f'Code Note Done yet')
    #     #print("flag2")
    # elif po is not None and pc is not None:
    #     pf = (pc*(Ae/C_star) - Ao*K_o*po) / (Ah*K_rp)
    #     pf_box.set_val(f'Code Not Done Yet')
    #     #print("flag3")                     # FIX!!!!!!!!!

    # plot the operating point
    p_ox_target = 437.5 # psi
    p_fu_target = 437.5 # psi
    if po is not None and pf is not None:
        if point_handle:
            point_handle.remove()
        main_ax.plot(po, pf, 'o', ms=10,
                                     color='white', mec='black', zorder=5)

        # compute OF ratio
        
        of_val = OF0[int(po/5), int(pf/5)]

        # T = m.Find_Temp(pc, of_val)
        # print("Temp:", T[1], ' K')

        # compute unchoke thresholds (Pc at which flow JUST unchokes)
        #	choked if P_up/Pc >= crit_ratio  →  unchokes when Pc > P_up/crit_ratio
        vapor_p_rp = 2000 / psi_to_pa #psi
        rp_p_loss = (0.5*(rho_rp) * 27.2**2)  / psi_to_pa   #rp-1 dynamic rpessure loss
        thresh_fuel = vapor_p_rp - rp_p_loss
        thresh_ox   = po / (max_stiff+1)
        thresh_pc   = min(thresh_fuel, thresh_ox)   # whichever orifice unchokes first
        # update text: OF and red-line Pc (psi)
        #out = [0,0,0,0,0,0]
        #print('CEA->', m.RunCEA(pc/psi_to_pa, 30, 'RP-1', 1, 298, 'O2', 298, of_val, 0, 0, 0, 1, 0, 'pythonCEA'))
        #print("out:", out)
        #print('shape:', out.shape)
        

        of_text.set_text(
            f'OF = {of_val:.2f}    '
            f'Red-line Pc = {thresh_ox:.1f} psi'
        )
        plt.draw()

def update(_):
    """Re-draw contours when any diameter changes."""
    global cont_pc, pc_labels, cont_of, of_labels, hat
    dh, dox, de = map(float, (tb_h.text, tb_o.text, tb_e.text))
    Pc, OF, mask, mask_error = compute_fields(dh, dox, de)

    pcm.set_array(Pc.ravel())
    pcm.set_clim(Pc.min(), Pc.max())
    cbar.update_normal(pcm)

    if hasattr(hat, 'collections'):
        for grp in (cont_pc.collections, cont_of.collections, hat.collections):
            for coll in grp:
                coll.remove()
    for txt in (*pc_labels, *of_labels):
        txt.remove()

    cont_pc   = main_ax.contour(X, Y, Pc, levels=pick_levels(Pc),
                                colors='k', linewidths=1)
    pc_labels = main_ax.clabel(cont_pc, fmt='%d psi', inline=True, fontsize=10)

    cont_of   = main_ax.contour(X, Y, OF, levels=[0.4, 0.6, 0.8, 1.0],
                                colors='k', linewidths=1)
    of_labels = main_ax.clabel(cont_of, fmt='OF: %.2f', inline=True, fontsize=10)

    hat = main_ax.contourf(X, Y, mask, levels=[0.5,1],
                           colors=['none'], hatches=['\\\\'], zorder=2)
    if hasattr(hat, 'collections'):
        for coll in hat.collections:
            coll.set(edgecolor='red', linewidth=0.5)

    hat2 = main_ax.contourf(X, Y, mask_error, levels=[0.5,1],
                            colors=['none'], hatches = ['\\\\'], zorder=3)

    plt.draw()

# hook up callbacks
for box in (tb_h, tb_o, tb_e):
    box.on_submit(update)
for box in (pf_box, po_box, pc_box):
    box.on_submit(compute_inlet)
    #print("flag5")


update(None)       # first draw
plt.show()