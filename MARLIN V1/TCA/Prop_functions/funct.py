from rocketcea.cea_obj import CEA_Obj
import matplotlib.pyplot as plt
import numpy as np

def get_isp_heatmap(Fuel,Ox, mr_range, p_range, pamb):    
    CEA = CEA_Obj(oxName=Ox, fuelName=Fuel)
    Mr = np.arange(mr_range[0], mr_range[1], 0.1)
    Pin = np.arange(p_range[0], p_range[1], 1)
    isp_grid = np.zeros((len(Mr), len(Pin)))

    for i in range(len(Mr)):
        for j in range(len(Pin)):
            mr = Mr[i]
            pc = Pin[j]
            Eps = CEA.get_eps_at_PcOvPe(Pc=pc, MR=mr, PcOvPe=(pc / pamb))
            isp = CEA.estimate_Ambient_Isp(Pc=pc, MR=mr, eps=Eps, Pamb=pamb)
            mr_stoich = CEA.getMRforER(ERphi=1.0)
            isp_grid[i, j] = isp[0]

    return Mr, Pin, isp_grid, mr_stoich

def get_tempstar_heatmap(Fuel,Ox, mr_range, p_range, pamb):  
    rankineToKelvin = 5.0 / 9.0  
    CEA = CEA_Obj(oxName=Ox, fuelName=Fuel)
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
            tempstar_grid[i, j] = tempstar[0] * rankineToKelvin

    return Mr, Pin, tempstar_grid, mr_stoich