from rocketcea.cea_obj_w_units import CEA_Obj
import numpy as np


def _mixture_density(mr, rho_fuel, rho_ox):
    return mr / (mr + 1.0) * rho_ox + 1.0 / (mr + 1.0) * rho_fuel

def get_isp_heatmap(Fuel, Ox, mr_range, p_range, pamb, rho_fuel, rho_ox):
    CEA = CEA_Obj(
        oxName=Ox,
        fuelName=Fuel,
        pressure_units="bar",
        temperature_units="K",
        cstar_units="m/s",
        sonic_velocity_units="m/s",
        density_units="kg/m^3",
        isp_units="sec",
    )
    Mr = np.arange(mr_range[0], mr_range[1], 0.1)
    Pin = np.arange(p_range[0], p_range[1], 1)
    isp_grid = np.zeros((len(Mr), len(Pin)))
    density = np.zeros(len(Mr))

    for i, mr in enumerate(Mr):
        density[i] = _mixture_density(mr, rho_fuel, rho_ox)
        for j, pc in enumerate(Pin):
            Eps = CEA.get_eps_at_PcOvPe(Pc=pc, MR=mr, PcOvPe=(pc / pamb))
            isp = CEA.estimate_Ambient_Isp(Pc=pc, MR=mr, eps=Eps, Pamb=pamb)
            mr_stoich = CEA.getMRforER(ERphi=1.0)
            isp_grid[i, j] = isp[0]

    return Mr, Pin, isp_grid, mr_stoich,density


def get_tempstar_heatmap(Fuel, Ox, mr_range, p_range, pamb):
    CEA = CEA_Obj(
        oxName=Ox,
        fuelName=Fuel,
        pressure_units="bar",
        temperature_units="K",
        cstar_units="m/s",
        sonic_velocity_units="m/s",
        density_units="kg/m^3",
        isp_units="sec",
    )
    Mr = np.arange(mr_range[0], mr_range[1], 0.1)
    Pin = np.arange(p_range[0], p_range[1], 1)
    tempstar_grid = np.zeros((len(Mr), len(Pin)))

    for i in range(len(Mr)):
        for j in range(len(Pin)):
            mr = Mr[i]
            pc = Pin[j]
            Eps = CEA.get_eps_at_PcOvPe(Pc=pc, MR=mr, PcOvPe=(pc / pamb))
            tempstar = CEA.get_Temperatures(Pc=pc, MR=mr, eps=Eps)
            mr_stoich = CEA.getMRforER(ERphi=1.0)
            tempstar_grid[i, j] = tempstar[1]

    return Mr, Pin, tempstar_grid, mr_stoich
