from sympy import symbols, log, diff, sin, cos
from sympy.parsing.sympy_parser import parse_expr

ERROR = 1E-6

x = symbols('x')

func = parse_expr(input('f(x): '))
x_0 = float(input('x_0: '))
dec = diff(func, x)

contador = 0
fx_1 = 100

while abs(fx_1) >= ERROR:
    fx_0 = func.subs(x, x_0).evalf()
    decx_0 = dec.subs(x, x_0).evalf()
    x_1 = (x_0) - (fx_0 / decx_0)
    fx_1 = func.subs(x, x_1).evalf()
    e = abs((x_1 - x_0) / x_1)
    print('iteración:', contador)
    print('x_0:', x_0)
    print('fx_0:', fx_0)
    print('x_1:',x_1)
    print('fx_1:', fx_1)
    print('*' * 15)
    if e > ERROR:
        x_0 = x_1
    contador += 1
print('Raíz: ', x_1)
print('f(Raíz):', fx_1)