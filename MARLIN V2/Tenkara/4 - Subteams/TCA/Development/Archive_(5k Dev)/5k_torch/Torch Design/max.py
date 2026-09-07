from rocketcea.cea_obj import CEA_Obj

# Constants
g = 9.8
mix_eff = 0.8
feetToMeters = 0.3048
p_amb = 14.7

# ________Setpoints_________

# Running Stoich GH2/GOx:
of = 8
p_c = 140
thrust = 15

# Function to translate of, pc, thrust to a throat diameter
def GH2GOx_throat_sizing_function(of, pc, thrust):
    # Create a CEA Rocket Problem
    cea = CEA_Obj(oxName = "GOX", fuelName = "GH2")

    # Get Cstar
    cstar = cea.get_Cstar(Pc=pc, MR=of) * feetToMeters

    # Use CEA to get exhaust velocity
    # Can't use get_Throat_Isp() because it returns vacuum isp, not ambient
    print(cea.estimate_Ambient_Isp(Pc=pc, MR=of, eps=1, Pamb=p_amb, frozen=1))
    CF = cea.get_PambCf(Pamb = p_amb, Pc=pc, MR = of, eps=1)[0]
    isp = cstar * CF / 9.81
    print(cea.get_full_cea_output(Pc= pc, MR= of, eps=1))

    print(f"ISP(s): {isp}")
    v_e = isp * mix_eff

    print(f"Exhaust Velocity: {v_e}")



def main():
    GH2GOx_throat_sizing_function(of, p_c, thrust)

if __name__ == "__main__":
    main()