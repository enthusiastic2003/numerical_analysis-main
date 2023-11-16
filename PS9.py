import toolbox.BVPSolver as BVPS
import toolbox.LinSolver as LS
import math as m
def p(x):
    return (-2/x)

def q(x):
    return (2/(x**2))

def r(x):
    return (m.sin(m.log(x))/(x**2))

def y(x):
    return ((1.139)*x)+(-0.039/(x**2))-(0.3*m.sin(m.log(x)))-(0.1*m.cos(m.log(x)))

print("Actual y(2)=", y(2))
(yp,yps)=BVPS.LinearShooting(p,q,r,1,2,1,2,0.1,10)
print("Predicted y(2)=",yp)

def ypp(X,Y,YP):
    return (1/8)*(32+(2*(X**3))-Y*YP)

def fy(X,Y,YP):
    return (-1/8)*YP

def fyp(X,Y,YP):
    return (-1/8)*Y

(xs,ys,ypps)=BVPS.shoot_nonlinear(1,3,17,(43/3),20,0.00001,100,ypp,fy,fyp)
print("Soln0=",ys)

