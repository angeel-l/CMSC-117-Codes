import methods as sp
import numpy as np

L = np.array([[1,0,0],[1,2,0],[1,1,1]])

U = L.transpose()
b = np.array([1,3,3])

# print("Lower:\n", L)
# print("\nUpper:\n", U)
# print("\nb:\n", b)

x = sp.ForwardSubRow(L,b)
print("x = ",x)
y = sp.ForwardSubCol(L,b)
print("y = ",y)
z = sp.BackwardSubRow(U,b)
print("z = ",z)
s = sp.BackwardSubCol(U,b)
print("s = ",s)