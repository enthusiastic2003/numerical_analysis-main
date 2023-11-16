import toolbox.BVPSolver as BVPS
import math as m

def p(x):
    return (-2/(x))

def q(x):
    return (2/(x**2))

def r(x):
    return (m.sin(m.log(x)))/(x**2)

def f(X,Y,YP):
    return ((-2/X)*YP)+((2/X**2)*Y)+((m.sin(m.log(X)))/(X**2))

def fy(X,Y,P):
    return (2/(X**2))

def fyp(X,Y,YP):
    return (-2/(X**2))

(xs,ys)=BVPS.FiniteDifferenceLinear(1,2,1,2,p,q,r,10)
(YSS,YSSS,YSSSS)=BVPS.shoot_nonlinear(1,2,1,2,9,0.0001,100,f,fy,fyp)
print("XS=",xs)
print("YS=",ys)
print("Actual YS=",YSSS)


