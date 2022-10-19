import numpy as np
import methods as md

# number of entries
n = 100

# creating the arrays a, b, and c
a = np.ones(n, dtype=float)
b = np.arange(1,n+1,dtype=float)
c = np.empty(n,float)
c[::2] = 1
c[1::2] = 3

# function to create the upper triangular matrix
def block(n, a, b, c):
    U = np.zeros((n,n), dtype=float)
    U = U + np.diag(a)
    U = U + np.diag(b[:n-1], 1)
    U = U + np.diag(c[:n-2], 2)
    U = U + np.diag(a[:n-3], 3)

    return U

# Creating text files
file1a = open("Item1a.text", "w+")
file1b = open("Item1b.text", "w+")

# Code for 1.a
A = block(100, a,b,c)
br, norm1a = md.BackwardSubRow(A, a)
bc, norm1b = md.BackwardSubCol(A, a)

# Writing to file Item1a.txt
answers = ["Backward Substitution by Rows: \n\n%s\n" % br, "\nResidual Max Norm: %s" %  norm1a,
           "\nBackward Substitution by Columns: \n\n%s\n" % bc, "\nResidual Max Norm: %s" %  norm1b]
file1a.writelines(answers)


# Code for 1.b
A = A.transpose()
fr, norm2a = md.ForwardSubRow(A, b)
fc, norm2b = md.ForwardSubCol(A, b)

# Writing to file Item1b.txt
answers = ["Forward Substitution by Rows: \n\n%s\n" % fr, "\nResidual Max Norm: %s" %  norm2a,
           "\nForward Substitution by Columns: \n\n%s\n" % fc, "\nResidual Max Norm: %s" %  norm2b]
file1b.writelines(answers)