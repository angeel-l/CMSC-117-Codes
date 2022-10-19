import numpy as np
import methods as md

# Creating text files
file2a = open("Item2a.text", "w+")
file2b = open("Item2b.text", "w+")

# Code for 2.a
A = np.matrix([[50, 107, 36],[25,54,20],[31,66,21]], dtype = float)
n = len(A)
b = np.ones(n)
x, normx1 = md.LUSolve(A.getA(), b)

# Writing to file Item2a.txt
answers = ["Solution of matrix A: %s\n" % x, "\nResidual Max Norm: %s" %  normx1]
file2a.writelines(answers)

# Code for 2.b
A = np.matrix([[10,2,1],[2,20,-2],[-2,3,10]], dtype = float)
x, normx2 = md.LUSolve(A.getA(), b)

# Writing to file Item2b.txt
answers = ["Solution of matrix A: %s\n" % x, "\nResidual Max Norm: %s" %  normx2]
file2b.writelines(answers)