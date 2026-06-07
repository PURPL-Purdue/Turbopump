import pandas as pd
import numpy as np

def gamma_lookup(data_Path, OF, reac_names):
    df = pd.read_excel(data_Path, sheet_name=f"{reac_names[0]} & {reac_names[1]}")
    of_to_gamma = dict(zip(df["OF"], df["GAMMAs_throat"]))
    return of_to_gamma[OF]
