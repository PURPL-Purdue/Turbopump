import pandas as pd
import matplotlib.pyplot as plt

def plot_all_variables(
    excel_file="properties.xlsx",
    x_column="x_m",
    highlight_points=None,
    figsize=(8, 5)
):
    """
    Plot every column in the Excel file against x_column.

    Parameters
    ----------
    excel_file : str
        Path to Excel file.
    x_column : str
        Column to use as x-axis.
    highlight_points : dict
        Dictionary of the form:
        {
            "column_name": [(x1, y1), (x2, y2), ...],
            ...
        }
    figsize : tuple
        Figure size.
    """

    df = pd.read_excel(excel_file)

    if x_column not in df.columns:
        raise ValueError(f"'{x_column}' not found in dataframe.")

    x = df[x_column]

    if highlight_points is None:
        highlight_points = {}

    for column in df.columns:

        # Skip x-axis column
        if column == x_column:
            continue

        y = df[column]

        plt.figure(figsize=figsize)
        plt.scatter(x, y, s=20, label="Data")

        # Plot highlighted points if provided
        if column in highlight_points:
            pts = highlight_points[column]

            hx = [p[0] for p in pts]
            hy = [p[1] for p in pts]

            plt.scatter(
                hx,
                hy,
                s=150,
                c="red",
                edgecolors="black",
                label="Known Points",
                zorder=10
            )

        plt.xlabel("Position [m]")
        plt.ylabel(column)
        plt.title(f"{column} vs {x_column}")
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()


# Example usage
highlight_points = {
    "T_chamber [K]": [
        (-0.124284824, 2829.57111),
        (0, 2591.87239),
        (0.121705411, 1428.01846)
    ],
    "P_chamber [bar]": [
        (-0.124284824, 40.02256),
        (0, 22.96619),
        (0.121705411, 0.82737)
    ],
    "gamma": [
        (-0.124284824, 1.18922),
        (0, 1.19985),
        (0.121705411, 1.22801)
    ],
    "Cp [(KJ/kg-K)]": [
        (-0.124284824, 2.85422226),
        (0, 2.60180862),
        (0.121705411, 2.20932387)
    ],   

}

plot_all_variables(
    excel_file="properties.xlsx",
    x_column="x_m",
    highlight_points=highlight_points
)