from sympy import *
from sympy.printing.cxxcode import cxxcode

x = symbols('x:3')
y = symbols('y:3')
a = symbols('a')
b = symbols('b')
c = symbols('c')


r_12 = Function('r_12')(x,y)
r_23 = Function('r_23')(x,y)
r_31 = Function('r_31')(x,y)

r_12 = sqrt((x[0]-x[1])**2 + (y[0]-y[1])**2)
r_23 = sqrt((x[1]-x[2])**2 + (y[1]-y[2])**2)
r_31 = sqrt((x[2]-x[0])**2 + (y[2]-y[0])**2)

V = Function('V')(r_12,r_23,r_31)

V = (add function to)

init_printing()

#Printing partial derivatives

#print(simplify(diff(V, x[0])))

for i in range(3):
    #print((diff(V, x[i])))
    #print((diff(V, y[i])))

    print(cxxcode(diff(V, x[i])))
    print(cxxcode(diff(V, y[i])))

    #print(latex((diff(V, x[i]))))
    #print(latex((diff(V, y[i]))))
