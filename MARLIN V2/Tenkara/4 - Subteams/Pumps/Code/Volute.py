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
    b3: Q_[float] # Volute width at impeller outlet
    h_gap: Q_[float] # Volute gap height at impeller outlet

lox_b3 = Q_(2.4, "millimeters")
ipa_b3 = Q_(2.4, "millimeters")
h_gap = Q_(1, "millimeters")
clearance = Q_(0.008, "inches")

def getKv(N_s: float) -> float:
    a, b = emp.Kv()
    K_v = a * N_s**b
    return K_v

def getVoluteVelocity(K_v: float, H: Q_[float]) -> Q_[float]:   # Uses constant velocity Stepanoff method
    avg_V = K_v * np.sqrt(2 * g.to('ft/s**2') * H.to('feet'))
    return avg_V

def getThroatArea(Q: Q_[float], avg_V: Q_[float]) -> Q_[float]:
    At = Q_(Q / avg_V, "in^2")
    return At

def getMinAngle(b3: Q_[float], h_gap: Q_[float], At: Q_[float]) -> Q_[float]:
    theta_min = (((b3 * h_gap) - ((b3**2) / 4) + ((b3**2) / 2)) / At) * Q_(360, "degrees")
    return theta_min

def getRhombusSideLength(At: Q_[float], theta: Q_[float], b3: Q_[float], h_gap: Q_[float]) -> Q_[float]:
    a = np.sqrt((At * (theta / Q_(360, "degrees"))) + ((b3**2) / 4) - (b3 * h_gap))
    return a

def getAreaAtAngle(Q: Q_[float], avg_V: Q_[float], e: Q_[float], e_sp: Q_[float], theta_min: Q_[float], b3: Q_[float], h_gap: Q_[float], n_volutes: int) -> Q_[float]:  # Get areas as a function of angle
    #if e < theta_min:
    #    A = (b3 * h_gap) + (e / theta_min) * ()
    #    return A
    #else:
    #    A = (e / e_sp) * (Q / (n_volutes * avg_V))
    A = (e / e_sp) * (Q / (n_volutes * avg_V))
    return A

def getVoluteTipHeight(a: Q_[float], b3: Q_[float], h_gap: Q_[float], impeller_r: Q_[float], clearance: Q_[float]) -> Q_[float]:
    h_top = (2 * a / np.sqrt(2)) - b3 / 2 + h_gap + impeller_r + clearance

    if h_top < 1.05 * impeller_r:
        print(f"!!! Below clearance !!!")

    return h_top

def calcRadialThrust(H: Q_[float], D2: Q_[float], B2_total: Q_[float], rho: Q_[float], Q: Q_[float], Q_bep: Q_[float]) -> Q_[float]:
    K = 0.36*(1 - (Q.to('gallons per minute') / Q_bep.to('gallons per minute'))**2)
    radial_thrust = (K * H * D2 * B2_total * rho * g).to('N')
    # [Charlie]: you don't have to convert the units each time, you can just wrap it all in parenthesis and convert once
    return radial_thrust

def sizeDiffuserLength(At: Q_[float], D_target: Q_[float], diffuserTheta: float) -> Q_[float]:
    Dt = np.sqrt(4 * At / math.pi)
    L = (D_target - Dt) / (2 * np.tan(diffuserTheta))
    return L

def calcDiffuserPressures(stagnationP1: Q_[float], At: Q_[float],  targetA: Q_[float], avg_V: Q_[float], rho: Q_[float]) -> Q_[float]:
    V2 = At / targetA * avg_V    
    staticP2 = stagnationP1 - (1/2 * rho * V2**2)
    return staticP2
    
class Volute:
    def __init__(self, geometry: VoluteGeometry, design_point: VoluteDesignPoint):
        self.K_v = getKv(self.N_s)
        self.avg_V = getVoluteVelocity(self.K_v, self.H)
        self.n_volutes = 2
        self.theta_min = getMinAngle(self.b3, self.h_gap, self.At)

        for theta in np.linspace(0, self.e_sp.magnitude, 100):
            self.A_theta = getAreaAtAngle(self.Q, self.avg_V, Q_(theta, 'radians'), self.e_sp, self.theta_min, self.b3, self.h_gap, self.n_volutes)

        self.diffuserLength = sizeDiffuserLength(self.At, self.D_target, self.diffuserTheta)
        self.staticP2 = calcDiffuserPressures(self.stagnationP1, self.At, self.targetA, self.avg_V, self.rho)

if __name__ == "__main__":
    lox_pump = Main.MakeDesign(Main.Propellant.LOX)
    lox_feed_P = lox_pump.DP.feed_pressure
    lox_H = lox_pump.DesignHeadRise
    lox_D2 = lox_pump.Impeller.d_2
    lox_B2_total = lox_pump.Impeller.b_2
    lox_rho = lox_pump.DP.rho
    lox_Q = lox_pump.DesignFlowRate
    lox_N_s = lox_pump.n_s
    lox_N_q = lox_pump.n_q
    _, lox_cm2, lox_cu2, _ = lox_pump.Impeller.GetOutletVelocities()
    lox_c2 = np.sqrt(lox_cm2**2 + lox_cu2**2)
    lox_K_v = getKv(lox_N_s)
    lox_avg_V = getVoluteVelocity(lox_K_v, lox_H)
    lox_At = getThroatArea(lox_Q, lox_avg_V)
    lox_theta_min = getMinAngle(lox_b3, h_gap, lox_At)
    lox_A_array = []
    lox_thetas = []
    for theta in np.linspace(0, 360, 361):
        A_theta = getAreaAtAngle(lox_Q, lox_avg_V, Q_(theta, 'degrees'), Q_(360, 'degrees'), lox_theta_min, lox_b3, h_gap, 1)
        lox_A_array.append(A_theta.to('cm^2').m)
        if A_theta == Q_(1, 'mm^2'):
            lox_thetas.append(theta)

    lox_stagnationP1 = lox_H * lox_rho * g
    lox_impellerStaticP = Q_(lox_stagnationP1 - (1/2 * lox_rho * lox_c2**2), "bar")
    lox_staticP1 = Q_(lox_stagnationP1 - (1/2 * lox_rho * lox_avg_V**2), "bar")
    lox_targetA = Q_(math.pi * ((12/16)/2)**2, "in^2")
    lox_throatA = Q_(lox_A_array[-1], "cm^2")
    lox_staticP2 = calcDiffuserPressures(lox_stagnationP1, lox_throatA, lox_targetA, lox_avg_V, lox_rho)
    lox_diffuserL = sizeDiffuserLength(lox_throatA, Q_(12/16, "inches"), Q_(7.5, "degrees")) # Placeholder diffuser length

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
    print(f"Nq metric specific speed: {lox_N_q:.4f}")
    print(f"LOX pump head rise: {lox_H.to('m'):.4f}")
    print(f"LOX pump diameter: {lox_D2.to('m'):.4f}")
    print(f"Lox pump width: {lox_B2_total.to('m'):.4f}")
    print(f"K_v: {lox_K_v:.4f}")
    print(f"Impeller static pressure: {lox_impellerStaticP.to('bar'):.4f}")
    print(f"Average velocity: {lox_avg_V.to('m/s'):.4f}")
    print(f"Tongue placement angle: {lox_theta_min.to("degrees"):.4f}")
    print(f"Start area: {Q_(lox_A_array[0], "cm^2").to('mm^2'):.4f}")
    print(f"Throat area: {Q_(lox_A_array[-1], "cm^2").to("in^2"):.4f}")
    print(f"Static pressure inside LOX volute: {lox_staticP1.to('bar'):.4f}")
    print(f"Static pressure at LOX diffuser exit: {lox_staticP2.to('bar'):.4f}")
    print(f"Diffuser length: {lox_diffuserL.to('millimeters'):.4f}")

    ipa_pump = Main.MakeDesign(Main.Propellant.IPA)
    ipa_feed_P = ipa_pump.DP.feed_pressure
    ipa_H = ipa_pump.DesignHeadRise
    ipa_D2 = ipa_pump.Impeller.d_2
    ipa_B2_total = ipa_pump.Impeller.b_2
    ipa_rho = ipa_pump.DP.rho
    ipa_Q = ipa_pump.DesignFlowRate
    ipa_N_s = ipa_pump.n_s
    ipa_N_q = ipa_pump.n_q
    _, ipa_cm2, ipa_cu2, _ = ipa_pump.Impeller.GetOutletVelocities()
    ipa_c2 = np.sqrt(ipa_cm2**2 + ipa_cu2**2)
    ipa_K_v = getKv(ipa_N_s)
    ipa_avg_V = getVoluteVelocity(ipa_K_v, ipa_H)
    ipa_At = getThroatArea(ipa_Q, ipa_avg_V)
    ipa_theta_min = getMinAngle(ipa_b3, h_gap, ipa_At)
    ipa_A_array = []
    ipa_thetas = []
    for theta in np.linspace(0, 360, 361):
        A_theta = getAreaAtAngle(ipa_Q, ipa_avg_V, Q_(theta, 'degrees'), Q_(360, 'degrees'), ipa_theta_min, ipa_b3, h_gap, 1)
        ipa_A_array.append(A_theta.to('cm^2').m)
        if A_theta == Q_(1, 'mm^2'):
            ipa_thetas.append(theta)

    ipa_stagnationP1 = ipa_H * ipa_rho * g
    ipa_impellerStaticP = Q_(ipa_stagnationP1 - (1/2 * ipa_rho * ipa_c2**2), "bar")
    ipa_staticP1 = Q_(ipa_stagnationP1 - (1/2 * ipa_rho * ipa_avg_V**2), "bar")
    ipa_targetA = Q_(math.pi * ((12/16)/2)**2, "in^2")
    ipa_throatA = Q_(ipa_A_array[-1], "cm^2")
    ipa_staticP2 = calcDiffuserPressures(ipa_stagnationP1, ipa_throatA, ipa_targetA, ipa_avg_V, ipa_rho)
    ipa_diffuserL = sizeDiffuserLength(ipa_throatA, Q_(12/16, "inches"), Q_(7.5, "degrees")) # Placeholder diffuser length
    
    #stagnationP1 = ipa_H.to("m") * ipa_rho.to("kg/m^3") * g.to("m/s^2")
    #staticP1 = stagnationP1 - 1/2 * ipa_rho.to("kg/m^3") * avg_V.to("m/s")**2
    #ipa_targetA = Q_(math.pi * ((12/16)/2)**2, "in^2")
    #throatA = Q_(A_array[-1], "m^2")
    #staticP2 = calcDiffuserPressures(stagnationP1, throatA, ipa_targetA.to("m^2"), avg_V.to("m/s"), ipa_rho.to("kg/m^3"))

    #plt.plot(np.linspace(0, 360, 361), A_array)
    #plt.xlabel = 'Volute wrap angle (degrees)'
    #plt.ylabel = 'Volute cross-sectional area (cm^2)'
    #plt.title = 'Volute cross-sectional area vs. wrap angle'
    #plt.grid(True)
    #plt.show()

    print(f"Specific speed: {ipa_N_s:.4f}")
    print(f"Nq metric specific speed: {ipa_N_q:.4f}")
    print(f"IPA pump head rise: {ipa_H.to('m'):.4f}")
    print(f"IPA pump diameter: {ipa_D2.to('m'):.4f}")
    print(f"IPA pump width: {ipa_B2_total.to('m'):.4f}")
    print(f"K_v: {ipa_K_v:.4f}")
    print(f"Impeller static pressure: {ipa_impellerStaticP.to('bar'):.4f}")
    print(f"Average velocity: {ipa_avg_V.to('m/s'):.4f}")
    print(f"Tongue placement angle: {ipa_theta_min.to("degrees"):.4f}")
    print(f"Start area: {Q_(ipa_A_array[0], "cm^2").to('mm^2'):.4f}")
    print(f"Throat area: {Q_(ipa_A_array[-1], "cm^2").to("in^2"):.4f}")
    print(f"Static pressure inside IPA volute: {ipa_staticP1.to('bar'):.4f}")
    print(f"Static pressure at IPA diffuser exit: {ipa_staticP2.to('bar'):.4f}")
    print(f"Diffuser length: {ipa_diffuserL.to('millimeters'):.4f}")

    #radial_thrust = calcRadialThrust(
    #    lox_pump.DesignHeadRise,
    #    lox_pump.Impeller.d_2,
    #    lox_pump.Impeller.b_2,
    #    lox_pump.DP.rho,
    #    0 * lox_pump.DesignFlowRate,
    #    lox_pump.DesignFlowRate
    #    )
    #print(f"Radial_thrust: {radial_thrust.to('N'):.4f}")

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
