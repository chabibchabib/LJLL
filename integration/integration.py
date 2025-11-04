import sympy as sp
from math import exp
# Définir les variables
x, y,z = sp.symbols('x y z')

# Définir la fonction
f = x**5 + y**2 +2*x*y+z**2

# Integrale 3-simplex
I = sp.integrate(f, (y, 0, 1-x-z), (x, 0, 1-z), (z, 0, 1))

print("Intégrale symbolique :", I)
print("Valeur numérique :", float(I.evalf()))

ff=21*x+x**2-3
II = sp.integrate(ff, (x, 0, 1))
print("Valeur numérique 1D :", float(II.evalf()))