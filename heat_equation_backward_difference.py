import numpy as np
from typing import Callable


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

if __name__ == "__main__":
    f: Callable[[float], float] = lambda x: np.sin(np.pi * x)
    l: float = 1.0
    T: float = 1.0
    alpha: float = 1.0
    m: int = 10
    N: int = 100
    solution = heat_equation_backward_difference(f, l, T, alpha, m, N)
    print("solution: w_{i, 50}", solution[:, 50])
