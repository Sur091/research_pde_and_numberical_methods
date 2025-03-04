import numpy as np
from typing import Callable
import matplotlib.pyplot as plt


def heat_equation_backward_difference(f: Callable[[float], float], l:float, T:float, alpha: float, m: int, N: int) -> np.ndarray:
    """
    Solve the heat equation using the backward difference method.

    INPUT:
    - l: End Point
    - T: Maximum Time
    - alpha: Constant
    - m: `>=` 3
    - N: `>=` 1
    - f: Boundary condition such that `f(x) = u(x, 0)`

    OUTPUT:
    - Approxiamtions `w_{i,j}` to `u(x_i, t_j)` for each `i = 1,..., m-1` and `j = 1,..., N`
    """
    # STEP 1
    h: float = l / m
    k: float = T / N
    lam: float = alpha * alpha * k / (h * h)

    # Solution array
    W = np.zeros(m)
    solution = np.zeros((m, N+1))

    # ls
    l_array = np.zeros(m)
    u_array = np.zeros(m)
    z_array = np.zeros(m)

    # STEP 2
    for i in range(1, m):
        W[i] = f(i * h)

    # STEP 3
    l_array[1]= 1.0 + 2.0 * lam;
    u_array[1]= -lam / l_array[1];

    # STEP 4
    for i in range(2, m-1):
        l_array[i] = 1.0 + 2.0 * lam + lam * u_array[i-1];
        u_array[i] = -lam / l_array[i];

    # STEP 5
    l_array[m-1] = 1.0 + 2.0 * lam + lam * u_array[m-2]

    # STEP 6
    for j in range(1, N+1):
        # STEP 7
        t: float = j * k;
        z_array[1] = W[1] / l_array[1]

        # STEP 8
        for i in range(2, m):
            z_array[i] = (W[i] + lam * z_array[i-1]) / l_array[i]
        # print("z_array", z_array)

        # STEP 9
        W[m-1] = z_array[m-1]

        # STEP 10
        for i in range(m-2, 0, -1):
            W[i] = z_array[i] - u_array[i] * W[i+1]

        # STEP 11
        solution[:, j] = W.copy()
    return solution


def main():
    f: Callable[[float], float] = lambda x: np.sin(np.pi * x)
    l: float = 1.0
    T: float = 0.5
    alpha: float = 1.0
    m: int = 10
    N: int = 100
    # The solution is m x N+1 matrix
    solution = heat_equation_backward_difference(f, l, T, alpha, m, N)
    
    # Create x and t arrays for plotting
    x = np.linspace(0, l, m)         # m points for the spatial domain
    t = np.linspace(0, T, N+1)         # N+1 time points

    # Create a meshgrid. Note: meshgrid returns arrays of shape (len(t), len(x))
    X, T_grid = np.meshgrid(x, t)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="3d")
    
    # Set a new view for the camera: adjust elevation and azimuth as desired
    ax.view_init(elev=10, azim=30)
    
    # Transpose the solution array to match the dimensions of X and T_grid
    surf = ax.plot_surface(X, T_grid, solution.T, cmap="viridis")
    ax.set_xlabel("x")
    ax.set_ylabel("t")
    ax.set_zlabel("u(x,t)")
    fig.colorbar(surf, shrink=0.5, aspect=5)
    plt.title("Heat Equation Solution Visualization")
    plt.savefig("heat_equation.png")
    # plt.show()
    



if __name__ == "__main__":
    main()
