#name: Erelle Boubli
#Id: 324460443
#name:roi putterman
#Id: 314919010

import math
import sympy

epsilon=0.0001
def Pol(x):
    return x**3-x-1
def Pol1(x):
    return x**3-math.cos(x)
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
    fxr = polinom(xr0)
    fdxr = dfx(xr0)
    diff=100
    xr=xr0
    i=0
    while diff>epsilon:
        i+=1
        xr1=FXR(fxr,fdxr,xr)
        fxr = polinom(xr1)
        fdxr = dfx(xr1)
        diff=xr-xr1
        xr=xr1
        print(xr,fxr,fdxr,diff)
    print("num of iterations is: ",i,"\nthe solution is: ",xr)
NewtonRaphson(Pol,Range1,epsilon)

def secant_method (polinom,Range,epsilon):
    print("secant_method:")
    x0 = Range[0]
    x1 = Range[1]
    fx=polinom(x0)
    diff=1
    print("xi=",x0, " xi+1=",x1, " f(xi)=",fx)
    while diff>epsilon:
        nextx=(x0*polinom(x1)-x1*polinom(x0))/(polinom(x1)-polinom(x0))
        x0=x1
        fx= polinom(x0)
        diff=x1-nextx
        x1=nextx
        print("xi=",x0, " xi+1=",x1, " f(xi)=",fx)
    print("solution=", x1)


secant_method(Pol1, Range2, epsilon)

