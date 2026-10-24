#   TODO: 
#   Integrate this with the impeller sizing (automatically wooooo)
#   integrate auto cross section dimension calcs

#   Assumed single volute
#   Diamond cross section?
#   Consider joined discharge double volute to reduce radial forces on rotor
#   Gulich p. 403:
#   If nevertheless a 6-bladed impeller (for example) is to be employed 
#   it is advisable to reduce the wrap angle of the inner volute εi to
#   165° or 170° so that two blades never pass the cutwaters simultaneously

from pint import Quantity as Q_
from dataclasses import dataclass
from matplotlib import pyplot as plt
import numpy as np
import math
import sys
from pathlib import Path
#sys.path.insert(0, str(Path(__file__).resolve().parents[2]))	# Pumps folder, so the Code package imports when run directly
import Main #TODO uhhhh talk to @Charlie Zhu
import EmpiricalRelations as emp
# Constants
g = Q_(9.81, 'm/s**2')

@dataclass
class VoluteDesignPoint:
    N_shaft: Q_[float] # Shaft speed (rpm)
    Q: Q_[float] # Pump volumetric flow rate
    H: Q_[float] # Pump head rise
    N_s: float # Pump specific speed

@dataclass
class VoluteGeometry:
    e_sp: Q_[float] # Wrap angle
    n_volutes: int # Number of volutes

def getKv(N_s: float) -> float:
    a, b = emp.Kv()
    K_v = a * N_s**b
    return K_v

def getVoluteVelocity(K_v: float, H: Q_[float]) -> Q_[float]:   # Uses constant velocity Stepanoff method
    avg_V = K_v * np.sqrt(2 * g.to('ft/s**2') * H.to('feet'))
    return avg_V

def getAreaAtAngle(Q: Q_[float], avg_V: Q_[float], e: Q_[float], e_sp: Q_[float], n_volutes: int) -> Q_[float]:  # Get areas as a function of angle
    #if e < Q_(1, 'degrees'):
    #    A = (Q_(1, 'degrees') / e_sp) * (Q / (n_volutes * avg_V))
    #    print(e, A)
    #    return A
    A = (e / e_sp) * (Q / (n_volutes * avg_V))

    min_gap = Q_(1, "mm")
    if np.sqrt(A) < min_gap:
        A = min_gap**2
    return A

def calcRadialThrust(H: Q_[float], D2: Q_[float], B2_total: Q_[float], rho: Q_[float], Q: Q_[float], Q_bep: Q_[float]) -> Q_[float]:
    K = 0.36*(1 - (Q.to('gallons per minute') / Q_bep.to('gallons per minute'))**2)
    radial_thrust = (K * H * D2 * B2_total * rho * g).to('N')
    # [Charlie]: you don't have to convert the units each time, you can just wrap it all in parenthesis and convert once
    return radial_thrust

def sizeDiffuserLength(At: Q_[float], D_target: Q_[float], diffuserTheta: float) -> Q_[float]:
    Dt = math.sqrt(4 * At / math.pi)
    L = (D_target - Dt) / (2 * math.tan(diffuserTheta))
    return L

def calcDiffuserPressures(stagnationP1: Q_[float], At: Q_[float],  targetA: Q_[float], avg_V: Q_[float], rho: Q_[float]) -> Q_[float]:
    V2 = At / targetA * avg_V    
    staticP2 = stagnationP1 - (1/2 * rho * V2**2)
    return staticP2
    
class Volute:
    def __init__(self, geometry: VoluteGeometry, design_point: VoluteDesignPoint):
        self.K_v = getKv(self.N_s)
        self.avg_V = getVoluteVelocity(self.K_v, self.H)

        for theta in np.linspace(0, self.e_sp.magnitude, 100):
            self.A_theta = getAreaAtAngle(self.Q, self.avg_V, Q_(theta, 'radians'), self.e_sp, self.n_volutes)

        self.diffuserLength = sizeDiffuserLength(self.At, self.D_target, self.diffuserTheta)
        self.staticP2 = calcDiffuserPressures(self.stagnationP1, self.At, self.targetA, self.avg_V, self.rho)

if __name__ == "__main__":
    lox_pump = Main.CreatePump('LOX')
    lox_H = lox_pump.DesignHeadRise
    lox_D2 = lox_pump.Impeller.d_2
    lox_B2_total = lox_pump.Impeller.b_2
    lox_rho = lox_pump.DP.rho
    lox_Q = lox_pump.DesignFlowRate
    lox_N_s = lox_pump.n_s
    K_v = getKv(lox_N_s)
    avg_V = getVoluteVelocity(K_v, lox_H)
    plt.figure()
    A_array = []
    thetas = []
    for theta in np.linspace(0, 360, 361):
        A_theta = getAreaAtAngle(lox_Q, avg_V, Q_(theta, 'degrees'), Q_(360, 'degrees'), 1)
        A_array.append(A_theta.to('cm^2').m)
        if A_theta == Q_(1, 'mm^2'):
            thetas.append(theta)

    stagnationP1 = lox_H * lox_rho * g
    staticP1 = Q_(stagnationP1 - (1/2 * lox_rho * avg_V**2), "bar")
    lox_targetA = Q_(math.pi * ((12/16)/2)**2, "in^2")
    throatA = Q_(A_array[-1], "cm^2")
    staticP2 = calcDiffuserPressures(stagnationP1, throatA, lox_targetA, avg_V, lox_rho)
    
    #stagnationP1 = lox_H.to("m") * lox_rho.to("kg/m^3") * g.to("m/s^2")
    #staticP1 = stagnationP1 - 1/2 * lox_rho.to("kg/m^3") * avg_V.to("m/s")**2
    #lox_targetA = Q_(math.pi * ((12/16)/2)**2, "in^2")
    #throatA = Q_(A_array[-1], "m^2")
    #staticP2 = calcDiffuserPressures(stagnationP1, throatA, lox_targetA.to("m^2"), avg_V.to("m/s"), lox_rho.to("kg/m^3"))

    #plt.plot(np.linspace(0, 360, 361), A_array)
    #plt.xlabel = 'Volute wrap angle (degrees)'
    #plt.ylabel = 'Volute cross-sectional area (cm^2)'
    #plt.title = 'Volute cross-sectional area vs. wrap angle'
    #plt.grid(True)
    #plt.show()

    print(f"Specific speed: {lox_N_s:.4f}")
    print(f"LOX pump head rise: {lox_H.to('m'):.4f}")
    print(f"LOX pump diameter: {lox_D2.to('m'):.4f}")
    print(f"Lox pump width: {lox_B2_total.to('m'):.4f}")
    print(f"K_v: {K_v:.4f}")
    print(f"Average velocity: {avg_V.to('m/s'):.4f}")
    print(f"Tongue placement angle: {Q_(thetas[-1], "degrees").to("degrees"):.4f}")
    print(f"Start area: {Q_(A_array[0], "cm^2").to('mm^2'):.4f}")
    print(f"Throat area: {Q_(A_array[-1], "cm^2").to("in^2"):.4f}")
    print(f"Static pressure inside LOX volute: {staticP1.to('bar'):.4f}")
    print(f"Static pressure at LOX diffuser exit: {staticP2.to('bar'):.4f}")

    radial_thrust = calcRadialThrust(
        lox_pump.DesignHeadRise,
        lox_pump.Impeller.d_2,
        lox_pump.Impeller.b_2,
        lox_pump.DP.rho,
        0 * lox_pump.DesignFlowRate,
        lox_pump.DesignFlowRate
        )
    print(f"Radial_thrust: {radial_thrust.to('N'):.4f}")

# System inputs
#Q_design = Q_(1.78, 'liters/second')
#H_design = Q_(294.8, 'meters')

# Impeller inputs
#C_u2 =  Q_(34.516, 'm/s')   # tangential velocity component
#r2 = Q_(1.1, 'inches')      # impeller outlet radius
#C_2 = Q_(68.2, 'm/s')       # Absolute velocity

# Volute inputs
#e_sp = Q_(2 * math.pi, 'radians')          # volute wrap angle (2pi for single volutes)
#K_v = 0.55                                 # experimental average velocity coefficient from h&h equation 6-70 pg 220
                                           # and from https://ijsea.com/archive/volume8/issue8/IJSEA08081016.pdf 

# Throat area calculations
#avg_V = K_v * np.sqrt(2 * g.to('ft/s**2') * H_design.to('feet'))    # use empirical relation from h&h pg 220
#A_throat = Q_design / avg_V
#throat_diameter = np.sqrt((4 / math.pi) * A_throat)                 # if it was a circle
#HH_throatV = Q_design.to('gallons per minute') / (A_throat.to('in**2') * 3.12) # flow velocity at throat outlet according to H&H pg 222

#print(f"Throat area: {A_throat.to('cm**2'):.4f}")
#print(f"Throat diameter: {throat_diameter.to('cm'):.4f}")
#print(f"Average velocity: {avg_V.to('m/s'):.4f}")
#print(f"Throat velocity: {HH_throatV.to('m/s'):.4f}")

# From gulich pg 420:
#X_sp = (Q_design / (math.pi * C_u2 * r2)) * ((e_sp / Q_(2 * math.pi, 'radians')).to('dimensionless'))
