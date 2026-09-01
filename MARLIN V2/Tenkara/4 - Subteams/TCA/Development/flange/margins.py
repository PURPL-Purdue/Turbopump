def margins(Safety_Factor, Material_Yield_Stress, Max_Expected_Stress):
    sigma_required = Max_Expected_Stress * Safety_Factor
    margin = Material_Yield_Stress/sigma_required - 1
    return margin * 100