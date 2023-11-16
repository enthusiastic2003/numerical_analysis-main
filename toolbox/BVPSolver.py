import sys
sys.path.append('./toolbox/')

import Deg2Solver as D2S
import LinSolver as LS

def LinearShooting(p,q,r,alpha,beta,a,b,h,n):
    """ 
    This functions solves BVP, and returns the function value and its derivative at all the mesh
    points , using the linear shooting and RK4 method.
    PROBLEM:
    Given: ypp(t,y,yp)=p(t)*yp(t)+q(t)*y(t)+r(t) ;
           y(a)=alpha ;
           y(b)=beta ;
           a+h*n=b
    """
    def f1(X,Y,YP):
        return (p(X)*YP)+(q(X)*Y)+r(X)
    
    def f2(X,Y,YP):
        return (p(X)*YP)+(q(X)*Y)
    
    (y1is,y1pis)=D2S.Deg2RK4(f1,a,alpha,0,h,n)
    (y2is,y2pis)=D2S.Deg2RK4(f2,a,0,1,h,n)
    yps=[]
    ys=[]
    for i in range(n+1):
        yps.append(y1pis[i] + (((beta - y1is[n]) / y2is[n]) * y2pis[i]))
    for i in range(n+1):
        ys.append(y1is[i] + (((beta - y1is[n]) / y2is[n]) * y2is[i]))
    return (ys,yps)

def shoot_nonlinear(a,b,alpha, beta, n, tol, M,f,fy,fyp):
    """
    Given: f(x,y,yp) ; a<=x<=b ; y(a)=alpha ; y(b)=beta 
    b=a+nh;
    Tolerance=tol
    M is max number of allowed iterations
    fy = df/dy
    fyp = df/ dyp
    """

    w1 = [0 for i in range(n+1)]
    w2 = [0 for i in range(n+1)]
    h = (b-a)/n
    k = 1
    TK = (beta - alpha)/(b - a)

    while k <= M:

        w1[0] = alpha
        w2[0] = TK
        u1    = 0
        u2    = 1

        for i in range(1,n+1):
            x = a + (i-1)*h     #step 5

            t = x + 0.5*(h)

            k11 = h*w2[i-1]     #step 6

            k12 = h*f(x,w1[i-1],w2[i-1])
            k21 = h*(w2[i-1] + (1/2)*k12)
            k22 = h*f(t, w1[i-1] + (1/2)*k11, w2[i-1] + (1/2)*k12)
            k31 = h*(w2[i-1] + (1/2)*k22)
            k32 = h*f(t, w1[i-1] + (1/2)*k21, w2[i-1] + (1/2)*k22)
            t   = x + h
            k41 = h*(w2[i-1]+k32)
            k42 = h*f(t, w1[i-1] + k31, w2[i-1] + k32)
            w1[i] = w1[i-1] + (k11 + 2*k21 + 2*k31 + k41)/6
            w2[i] = w2[i-1] + (k12 + 2*k22 + 2*k32 + k42)/6   
            kp11 = h*u2
            kp12 = h*(fy(x,w1[i-1],w2[i-1])*u1 + fyp(x,w1[i-1], w2[i-1])*u2)
            t    = x + 0.5*(h)
            kp21 = h*(u2 + (1/2)*kp12)
            kp22 = h*((fy(t, w1[i-1],w2[i-1])*(u1 + (1/2)*kp11)) + fyp(x+h/2, w1[i-1],w2[i-1])*(u2 +(1/2)*kp12))
            kp31 = h*(u2 + (1/2)*kp22)
            kp32 = h*((fy(t, w1[i-1],w2[i-1])*(u1 + (1/2)*kp21)) + fyp(x+h/2, w1[i-1],w2[i-1])*(u2 +(1/2)*kp22))
            t    = x + h
            kp41 = h*(u2 + kp32)
            kp42 = h*(fy(t, w1[i-1], w2[i-1])*(u1+kp31) + fyp(x + h, w1[i-1], w2[i-1])*(u2 + kp32))
            u1 = u1 + (1/6)*(kp11 + 2*kp21 + 2*kp31 + kp41)
            u2 = u2 + (1/6)*(kp12 + 2*kp22 + 2*kp32 + kp42)


        r = abs(w1[n] - beta)
        #print(r)
        if r < tol:
                xs = [a + i*h for i in range(n+1)]
                return (xs,w1,w2)            

        TK = TK -(w1[n]-beta)/u1

        k = k+1


    print("Maximum number of iterations exceeded")   
    return
    
def FiniteDifferenceLinear(a,b,alpha,beta,p,q,r,n):
    """
    Applies Finite Differentiation Formulation to calculate the values of y at mesh points. 
    Creates the Matrices A and b, and then passes them to any of the Linear equation Solvers.
    """
    h=(b-a)/(n)

    xs=[a+i*h for i in range(0,n+1)]

    b=[(-1*(h**2)*r(xs[1]))+((1+h*0.5*p(xs[1]))*alpha)]
    for i in range(2,n-1):
        b.append((-1*(h**2)*r(xs[i])))
    b.append((-1*(h**2)*r(xs[n-1]))+((1-(h*p(xs[n-1])*0.5))*beta))
    A=[[0 for i in range(n-1)] for j in range(n-1)]

    A[0][0]=2+((h**2)*q(xs[1]))
    A[0][1]= -1+(h*0.5*p(xs[1]))

    A[n-2][n-2]=2+((h**2)*q(xs[n-1]))
    A[n-2][n-3]=-1-(h*0.5*p(xs[n-1]))

    for i in range(1,n-2):
        A[i][i]=2+((h**2)*q(xs[i+1]))
        A[i][i+1]=-1+(h*0.5*p(xs[i+1]))
        A[i][i-1]=-1-(h*0.5*p(xs[i+1]))
    w=LS.GaussElem(A,b)
    print("B=",b)
    return (xs,w)