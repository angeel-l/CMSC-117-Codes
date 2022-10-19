import numpy as np

def ForwardSubRow(L,b):
    """
    Solving Lx = b using forward substitution by rows.
    
    Parameters:
    L : matrix
        n x n lower triangular matrix
    b : vector
        n x 1   
        
    Output
    ------
    x : vector
        n x 1, solution to Lx = b 
    """
    n = len(L)
    x = np.zeros(n)
    x[0] = b[0]/L[0][0]
    for i in range(1,n):
        s = 0
        for j in range(0,i):
            s = s + L[i][j]*x[j]
        x[i] = (b[i]-s)/L[i][i]
    norm = MaxNorm(L, x, b)
    return (x, norm)

def ForwardSubCol(L,b):
    """
    Solving Lx = b using forward substitution by columns.
    
    Parameters:
    L : matrix
        n x n lower triangular matrix
    b : vector
        n x 1   
        
    Output
    ------
    x : vector
        n x 1, solution to Lx = b 
    """
    d = b
    n = len(L)
    for j in range(0,n-1):
        b[j] = b[j]/L[j][j]
        for i in range(j+1,n):
            b[i] = b[i] - L[i][j]*b[j]
    b[n-1] = b[n-1]/L[n-1][n-1]
    norm = MaxNorm(L, b, d)
    return (b, norm)

def BackwardSubRow(U,b):
    """
    Solving Ux = b using backward substitution by rows.
    
    Parameters:
    U : matrix
        n x n upper triangular matrix
    b : vector
        n x 1   
        
    Output
    ------
    x : vector
        n x 1, solution to Ux = b 
    """
    n = len(U)
    x = np.zeros(n)
    x[n-1] = b[n-1]/U[n-1][n-1] 
    for i in range(n-2,-1,-1):
        s = 0
        for j in range(i+1,n):
            s = s + U[i][j]*x[j]
        x[i] = (b[i]-s)/U[i][i]   
    norm = MaxNorm(U, x, b)
    return (x, norm) 
    
def BackwardSubCol(U, b):
    """
    Solving Ux = b using backward substitution by columns.
    
    Parameters
    ----------
    U : matrix
        n x n upper triangular matrix
    b : vector
        n x 1   
        
    Output
    ------
    b : Vector
        Solution to Ux = b.
    """
    d = b
    n = len(b)
    for j in range(n-1,0,-1):
        b[j] = b[j]/U[j][j]
        for i in range(0,j):
            b[i] -= U[i][j]*b[j]
    b[0] = b[0]/U[0][0]
    norm = MaxNorm(U, b, d)
    return (b, norm)

def MaxNorm(A, xtild, d):
    """
    Computing the Residual Max Norm.
    
    Parameters
    ----------
    A : matrix
        n x n triangular matrix
    xtild : vector
        n x 1 computed solution
    d : vector
        n x 1
    """
    sol = d - np.dot(A,xtild)
    rmn = 0
    for i in sol:
        rmn = rmn + abs(i)
    return rmn

def LUIJK(A):
    """
    Solving for the LU factorization in A stored in A

    Parameters:
    A: matrix
       n x n

    Output
    ------
    A: matrix
       n x n, LU factorization of the input A
    """
    n = len(A)
    for j in range(1,n):
        A[j][0] = A[j][0] / A[0][0]
    for i in range(1,n):
        for j in range(i,n):
            s = 0
            for k in range(0,i):
                s = s + A[i][k]*A[k][j]
            A[i][j] = A[i][j] - s
        for j in range(i+1, n):
            s = 0
            for k in range(0,i):
                s = s + A[j][k]*A[k][i]
            A[j][i] = (A[j][i] - s) / A[i][i]
    return A

def GetLU(A):
    """
    Solving for the LU factors of A

    Parameters:
    A: matrix
       n x n

    Output
    ------
    L: n x n, lower triangular matrix
    U: n x n, upper triangular matrix
    """
    A = LUIJK(A)
    n = len(A)
    L = np.zeros((n,n))
    U = np.zeros((n,n))
    for i in range(0,n):
        L[i][i] = 1
        for j in range(0,i):
            L[i][j] = A[i][j]
        for j in range(i,n):
            U[i][j] = A[i][j]
    return (L, U)
            
def LUSolve(A, b):
    """
    Solving for x by Ly = b and Ux = y
    
    Parameters:
    A: matrix
       n x n
    b: vector
       n x 1
    """
    L, U = GetLU(A)
    y, normf = ForwardSubRow(L, b)
    x, normb = BackwardSubRow(U, y)
    normx = MaxNorm(A, x, b)
    return x, normx
