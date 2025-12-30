# == PLOTS == 

def ring_points(radius_mm, N):
    ang = np.linspace(0, 2*np.pi, N, endpoint=False)
    x = radius_mm * np.cos(ang)
    y = radius_mm * np.sin(ang)
    return x, y

def plot_injector_layout(D_c_mm, R_ring_f_mm, R_ring_ox_mm, N_hole_ring, margin_extra_mm=5):
    combRad = D_c_mm / 2.0 #chamber radius
    x_f_in,  y_f_in  = ring_points(R_ring_f_mm[0], N_hole_ring)
    x_f_out, y_f_out = ring_points(R_ring_f_mm[1], N_hole_ring)
    x_ox_in,  y_ox_in  = ring_points(R_ring_ox_mm[0], N_hole_ring)
    x_ox_out, y_ox_out = ring_points(R_ring_ox_mm[1], N_hole_ring)
    fig, ax = plt.subplots(figsize=(10, 10))
    chamber = plt.Circle((0, 0), combRad, color='gray',
                         fill=False, linewidth=2)
    ax.add_artist(chamber)
    for R in R_ring_f_mm:
        ax.add_artist(plt.Circle((0, 0), R, color='red',
                                 fill=False, linewidth=0.7, alpha=0.4))
    for R in R_ring_ox_mm:
        ax.add_artist(plt.Circle((0, 0), R, color='blue',
                                 fill=False, linewidth=0.7, alpha=0.4))
    ax.scatter(x_f_in,  y_f_in, s=18, color='red', label='RP-1 inner')
    ax.scatter(x_f_out, y_f_out, s=18, color='darkred', label='RP-1 outer')
    ax.scatter(x_ox_in,  y_ox_in, s=18, color='royalblue', label='LOX inner')
    ax.scatter(x_ox_out, y_ox_out, s=18, color='navy', label='LOX outer')
    ax.set_aspect('equal', 'box')
    ax.set_xlabel('x [mm]')
    ax.set_ylabel('y [mm]')
    ax.set_title('4-ring impinging injector layout')
    ax.legend(loc='upper right')
    ax.grid(True, linestyle=':', linewidth=0.5)
    margin = combRad + margin_extra_mm
    ax.set_xlim(-margin, margin)
    ax.set_ylim(-margin, margin)
    plt.tight_layout()
    plt.show()

holes_per_ring = int(num_holes_lox_inj/Nrings)
plot_injector_layout(2*combRad * 1e3, Rring_rp1 * 1e3, Rring_lox * 1e3, holes_per_ring)