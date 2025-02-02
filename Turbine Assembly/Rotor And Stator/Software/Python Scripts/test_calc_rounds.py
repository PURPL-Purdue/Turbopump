import numpy as np
from scipy.optimize import fsolve

def solve_for_x1(c1, c2, r, y_2, x_2):
    """
    solves the equation (x+b)(r^2-(x+b)^2) = (y_2 - sqrt(r^2-(x+b)^2)) / (x_2 - x) for x.
    """
    initial_guess = - c1 - r
    def equation(x):
        left = - (x + c1) / np.sqrt(r**2 - (x + c1)**2)
        right = (y_2 - (c2 + np.sqrt(r**2 - (x + c1)**2))) / (x_2 - x)
        return left - right

    solution = fsolve(equation, initial_guess)
    
    y1 = np.sqrt(r**2 - (solution[0] + c1)**2) + c2

    return (solution[0], y1)

c1 = 0.00401292
c2 = 0.0026337
r = 0.0048
x2 = 0.06
y2 = 0.47566

print(solve_for_x1(c1, c2, r, y2, x2))