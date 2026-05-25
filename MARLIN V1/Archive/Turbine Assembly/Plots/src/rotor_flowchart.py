from diagrams import Diagram, Cluster, Edge
from diagrams.generic.blank import Blank

with Diagram("Turbine Design Flowchart", direction="LR"):

    with Cluster("Define Constants and Parameters"):
        A1 = Blank("Set Inlet Total Pressure, Temperature, Gas Properties")
        A2 = Blank("Initialize Mass Flow and Turbine RPM")
        A3 = Blank("Define Chord Length, Blade Spacing, Hub, and Rotor Radius")

    with Cluster("Calculate Turbine Geometry"):
        B1 = Blank("Calculate Throat Properties")
        B2 = Blank("calc_P_throat - Throat Pressure")
        B3 = Blank("calc_rho_throat - Throat Density")
        B4 = Blank("calc_v_throat - Throat Velocity")
        B5 = Blank("calc_A_throat - Throat Area")

    with Cluster("Back-Calculate Rotor Velocities"):
        C1 = Blank("Use rotorBackCalculate for Velocities, Angles")
        C2 = Blank("Calculate Relative and Absolute Velocities: v1, v2, w")
        C3 = Blank("Determine Flow Angles: a1, a2, b")
        C4 = Blank("plot_velocity_triangles_angles - Visualize Velocity Triangles")

    with Cluster("Generate Blade Geometry"):
        D1 = Blank("calculate_blade_efficiency - Blade Efficiency")
        D2 = Blank("generate_blade_geom - Define x/y Blade Profiles")
        D3 = Blank("Obtain Cross-Sectional Area, Max Thickness")
        D4 = Blank("plot_turbine - Plot Blade Geometry")

    with Cluster("Calculate Mach and Pressure Distribution"):
        E1 = Blank("calculateMachPressureDistribution - Mach/Pressure")
        E2 = Blank("Distribute Along Blade Surface")
        E3 = Blank("plotMachPressureDistributions - Visualize Distributions")

    with Cluster("Plot Turbine Geometry and Distributions"):
        F1 = Blank("Generate Complete Turbine Visualization")
        F2 = Blank("Overlay Mach/Pressure Results")

    # Connecting steps within clusters
    A1 >> A2 >> A3
    B1 >> B2 >> B3 >> B4 >> B5
    C1 >> C2 >> C3 >> C4
    D1 >> D2 >> D3 >> D4
    E1 >> E2 >> E3
    F1 >> F2

    # Connecting clusters horizontally
    A3 >> Edge(label="Next") >> B1
    B5 >> Edge(label="Next") >> C1
    C4 >> Edge(label="Next") >> D1
    D4 >> Edge(label="Next") >> E1
    E3 >> Edge(label="Next") >> F1
