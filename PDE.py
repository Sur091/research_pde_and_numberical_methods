import numpy as np


def finite_difference_poisson(f,g,a: int,b,c,d,n,m,TOL,iter):
    # Frequently used constants
    N = n + 1
    M = m + 1
    
    # Define mesh sizes
    h = (b - a)/n
    k = (d - c)/m

    # Define the mesh
    x = np.zeros(N)
    y = np.zeros(M)

    for i in range(N):
        x[i] = a + i*h

    for j in range(M):
        y[j] = c + j*k

    # Define the array of Function Values
    W = np.zeros((n-1, m-1))

    # Needed constants
    lam = (h/k)**2
    mu = 2*(1 + lam)
    ell = 1

    # Step 6 Gauss-Seidel iterations
    while ell <= iter:
        # Step 7
        z = (-(h**2)*f(x[1], y[m-1]) + g(x[0], y[m-1]) + lam*g(x[1], y[m]) + lam*W[0, m-3] + W[1, m-2])/mu
        norm = np.abs(z - W[0, m-2])
        W[0, m-2] = z

        # Step 8
        for i in range(2, n-1, 1):
            z = (-(h**2)*f(x[i], y[m-1]) + lam*g(x[i], y[m]) + W[i-2, m-2] + W[i, m-2] + lam*W[i-1, m-3])/mu
            # print("Step 8: z = ", z)

            if np.abs(W[i-1, m-2] - z) > norm:
                norm = np.abs(W[i-1, m-2] - z)
            
            W[i-1, m-2] = z

        

        # Step 9
        z = (-(h**2)*f(x[n-1], y[m-1]) + g(x[n], y[m-1]) + lam*g(x[n-1], y[m]) + W[n-3, m-2] + lam*W[n-2, m-3])/mu
        
        if np.abs(W[n-2, m-2] - z) > norm:
            norm = np.abs(W[n-2, m-2] - z)
        
        W[n-2, m-2] = z

        # Step 10
        for j in range(m-2, 1, -1):
            # Step 11
            z = (-(h**2)*f(x[1], y[j]) + g(x[0], y[j]) + lam*W[0, j] + lam*W[0, j-2] + W[1, j-1])/mu
            
            if np.abs(W[0, j-1] - z) > norm:
                norm = np.abs(W[0, j-1] - z)
            
            W[0, j-1] = z

            # Step 12
            for i in range(2, n-1, 1):
                z = (-(h**2)*f(x[i], y[j]) + W[i-2, j-1] + lam*W[i-1, j] + W[i, j-1] + lam*W[i-1, j-2])/mu

                if np.abs(W[i-1,j-1] - z) > norm:
                    norm = np.abs(W[i-1,j-1] - z)
                
                W[i-1,j-1] = z

            # Step 13
            z = (-(h**2)* f(x[n-1], y[j])+g(b, y[j]) + W[n-3, j-1] + lam * W[n-2, j] + lam* W[n-2, j-2]) / mu

            if np.abs(W[n-2, j-1] - z) > norm:
                norm = np.abs(W[n-2, j-1] - z)
            # norm = max(norm,) Change later
            
            W[n-2, j-1] = z
        
        # Step 14
        z = (-h**2 * f(x[1], y[1]) + g(a, y[1]) + lam * g(x[1], c) + lam * W[0, 1] + W[1, 0]) /mu

        if np.abs(W[0, 0] - z) > norm:
            norm = np.abs(W[0, 0] - z)
        W[0, 0] = z 

        # Step 15
        for i in range(2, n-1, 1):
            z = ( -h**2 * f(x[i], y[1]) + lam * g(x[i], c) + W[i-2, 0] + lam * W[i-1, 1] + W[i, 0] ) /mu

            if np.abs(W[i-1, 0] -z) > norm:
                norm = np.abs(W[i-1, 0] -z)
            W[i-1, 0] = z

        # Step 16
        z = (-h**2 * f(x[n-1], y[1]) + g(b, y[1]) + lam * g(x[n-1], c) + W[n-3, 0] + lam * W[n-2, 1]) / mu 
        if np.abs(W[n-2, 0] - z) > norm:
            norm = np.abs(W[n-2, 0] - z)
        W[n-2, 0] = z 

        # Step 17 onwards
        if norm <= TOL:
            return x, y, W # The procedure was successful
        else:
            ell += 1 
    print(W)
    print("Maximum number of iterations reached") # The procedure was unsuccessful


# efficiency: define h^2 as a variable
# efficiency: before each if statement define update = ... to avoid computing same thing twiceimport numpy as np

f = lambda x, y: 4

# Check for floating point equality
def g(x, y):
    if y == 0: return x**2
    if y == 2: return (x-2)**2
    if x == 0: return y ** 2
    if x == 1: return (y - 1) **2


a, b, c, d = 0, 1, 0, 2
n, m = 20, 40

TOL, iter = 1e-6, 1000

print(finite_difference_poisson(f, g, a, b, c, d, n, m, TOL, iter))
