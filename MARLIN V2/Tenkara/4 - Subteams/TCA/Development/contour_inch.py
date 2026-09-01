import pandas as pd
import numpy as np

MM2IN = 0.0393701

def main():
    file = "/Users/louis/Desktop/Turbopump_GIT/MARLIN V2/Tenkara/4 - Subteams/TCA/Outputs/contour.csv"    
    outfile = "/Users/louis/Desktop/Turbopump_GIT/MARLIN V2/Tenkara/4 - Subteams/TCA/Outputs\contour_inch.csv"
    df = pd.read_csv(file)
    x_in = df["x_mm"] * MM2IN
    y_in = df["y_mm"] * MM2IN

    df_in = pd.
    return



if __name__ == "__main__":
    main()