#name: Erelle Boubli
#Id: 324460443
#name:roi putterman
#Id: 314919010

import math
import sympy
import decimal

epsilon=0.0001
def Pol(x):
    return x**3-x-1
Range1=(1,2)
Range2=(0,1)

def FXR(fxr,fdxr,xr):
    return xr-(fxr/fdxr)

def dr(f):
    x = sympy.symbols('x')
    drx = sympy.diff(f)
    # dri=derivative(f,i)
    return drx

def dfx(i):
    x = sympy.symbols('x')
    f2 = Pol(x)
    f = dr(f2)
    f_num = f.subs([(x, i)])
    return f_num

def NewtonRaphson(polinom,Range,epsilon):
    print("NewtonRaphson:")
    a=Range[0]
    b=Range[1]
    xr0=(a+b)/2
    fxr = Pol(xr0)
    fdxr = dfx(xr0)
    diff=100
    xr=xr0
    i=0
    while diff>epsilon:
        i+=1
        xr1=FXR(fxr,fdxr,xr)
        fxr = Pol(xr1)
        fdxr = dfx(xr1)
        diff=xr-xr1
        xr=xr1
        print(xr,fxr,fdxr,diff)
    print("num of iterations is: ",i,"\nthe solution is: ",xr)
NewtonRaphson(Pol,Range1,epsilon)

def secant_method (polinom,Range,epsilon):
    print("secant_method:")
    a = Range[0]
    b = Range[1]


secant_method(Pol, Range2, epsilon)

