import pandas as pd
import matplotlib.pyplot as plt
import sys
import os


def plot_contours(*csv_files):
    if not csv_files:
        print("No CSV files provided.")
        return

    fig, ax = plt.subplots(figsize=(12, 5))

    for path in csv_files:
        df = pd.read_csv(path)
        label = os.path.splitext(os.path.basename(path))[0]

        # Sort by x so the line is drawn in order
        df = df.sort_values("x_mm")

        ax.plot(df["x_mm"], df["y_mm"], label=label, linewidth=1.5)

    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    ax.set_title("Nozzle Contour(s)")
    ax.legend()
    ax.set_aspect("equal")
    ax.grid(True, linestyle="--", alpha=0.5)
    plt.tight_layout()
    plt.show()


# python plot_csv.py file1.csv file2.csv ...
if __name__ == "__main__":
    files = sys.argv[1:]
    if not files:
        print("Usage: python plot_csv.py <file1.csv> [file2.csv ...]")
        sys.exit(1)
    plot_contours(*files)