import sympy
from sympy.matrices import Matrix

# I understand that this is possible: 
t = sympy.Symbol('t')
P0 = sympy.Matrix([sympy.Symbol('x1'), sympy.Symbol('y1')])
P1 = sympy.Matrix([sympy.Symbol('x2'), sympy.Symbol('y2')])
P2 = sympy.Matrix([sympy.Symbol('x3'), sympy.Symbol('y3')])
P3 = sympy.Matrix([sympy.Symbol('x4'), sympy.Symbol('y4')])

X = (1-t)**3*P0 + 3*t*(1-t)**2*P1 + 3*t**2*(1-t)*P2 + t**3*P3

Lx1 = sympy.Symbol('Lx1')
Lx2 = sympy.Symbol('Lx2')
Ly1 = sympy.Symbol('Ly1')
Ly2 = sympy.Symbol('Ly2')

print(X)

eqn1 = X[1]-Ly1-((X[0]-Lx1)/(Lx2-Lx1))*(Ly2-Ly1)
eqn2 = X[0]-Lx1-((X[1]-Ly1)/(Ly2-Ly1))*(Lx2-Lx1)

coeffs = sympy.Poly(eqn1, t).coeffs()
print(coeffs)

coeffs = sympy.Poly(eqn2, t).coeffs()
print(coeffs)