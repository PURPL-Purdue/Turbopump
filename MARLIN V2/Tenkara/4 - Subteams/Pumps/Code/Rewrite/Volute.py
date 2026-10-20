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
    N_shaft: Q_[float]
    Q: Q_[float]
    H: Q_[float]

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
    A = (e / e_sp) * (Q / (n_volutes * avg_V))
    return A

def calcRadialThrust(H: Q_[float], D2: Q_[float], B2_total: Q_[float], rho: Q_[float], Q: Q_[float], Q_bep: Q_[float]) -> Q_[float]:
    K = 0.36*(1 - (Q.to('gallons per minute') / Q_bep.to('gallons per minute'))**2)
    radial_thrust = K * H.to('ft') * D2.to('inches') * B2_total.to('inches') * rho.to('lb/ft**3') * g.to('ft/s**2')
    # [Charlie]: you don't have to convert the units each time, you can just wrap it all in parenthesis and convert once
    return radial_thrust

lox_pump = Main.CreatePump('LOX')
radial_thrust = calcRadialThrust(
    lox_pump.DesignHeadRise,
    lox_pump.Impeller.d_2,
    lox_pump.Impeller.b_2,
    lox_pump.DP.rho,
    lox_pump.DesignFlowRate,
    lox_pump.DesignFlowRate
    )

class Volute:
    def __init__(self, geometry: VoluteGeometry, design_point: VoluteDesignPoint):
        pass

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
