import numpy as np
from pint import Quantity as Q_
import plotly.graph_objects as go
from CoolProp import CoolProp as CP
from Pump import Pump, Propellant
import Pump1D
import Main

water_rho = Q_(CP.PropsSI('D', 'P', 101325, 'T', 298, 'Water'), 'kg/m**3')

def getWaterPower(Pump:Pump, N: Q_[float]) -> Q_[float]:
    prop_rho = Pump.DP.rho
    Q = (Pump.DP.mass_flow/Pump.DP.rho) * (N / Pump.DP.N_shaft)
    prop_power = Pump1D.GetPowerConsumption(Pump, Q, N)
    water_power = prop_power * (water_rho / prop_rho).to('dimensionless')
    return water_power

def getTorque(N: Q_[float], power: Q_[float]) -> Q_[float]:
    torque = (power / N).to('newton * meter')

    return torque

def sweepPowers(Pump:Pump, N_min: Q_[float], N_max: Q_[float]) -> Q_[np.ndarray]:
    powers = [getWaterPower(Pump, N) for N in np.linspace(N_min, N_max, 100)]
    
    return Q_(np.array([power.to('kW').m for power in powers]), 'kW')

def sweepTorque(N_array: Q_[np.ndarray], Power_array: Q_[np.ndarray]) -> Q_[np.ndarray]:
    return (Power_array / N_array).to('newton * meter')


if __name__ == "__main__":
    lox_pump = Main.CreatePump(Main.Propellant.LOX)
    ipa_pump = Main.CreatePump(Main.Propellant.IPA)
    N_min = Q_(1000, 'rpm')
    N_max = Q_(25000, 'rpm')
    lox_power_array = sweepPowers(lox_pump, N_min, N_max)
    ipa_power_array = sweepPowers(ipa_pump, N_min, N_max)
    speed_array = Q_(np.linspace(N_min.to('rpm').m, N_max.to('rpm').m, len(lox_power_array)),'rpm')
    lox_torque_array = sweepTorque(speed_array,lox_power_array)
    ipa_torque_array = sweepTorque(speed_array,ipa_power_array)

    fig = go.Figure(go.Scatter(
        x=speed_array,
        y=lox_power_array.to('kW').m,
        mode='lines',
        line=dict(color='blue', width=4),
        name='LOx Power Required'
    ))
    fig.update_layout(
        title='LOx Power Requirement vs. Shaft Speed',
        xaxis_title='Shaft speed [rpm]',
        yaxis_title='Power [kW]',
        template='plotly_white'
    )
    fig.show()

    fig = go.Figure(go.Scatter(
        x=speed_array,
        y=ipa_power_array.to('kW').m,
        mode='lines',
        line=dict(color='orange', width=4),
        name='IPA Power Required'
    ))
    fig.update_layout(
        title='IPA Power Requirement vs. Shaft Speed',
        xaxis_title='Shaft speed [rpm]',
        yaxis_title='Power [kW]',
        template='plotly_white'
    )
    fig.show()
