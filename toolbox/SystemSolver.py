import sys
sys.path.append('./toolbox/')

def System2eq_RK4(alpha1,alpha2,f1,f2,h,n,a):
    """ 
    This functions solves an IVP system of 2 equations and returns the value at all the mesh
    points of both the equations, with Runge-Kutta 4.
    PROBLEM:
    Given: f1(t,u1,u2)=u1'
           f2(t,u1,u2)=u2'
           u1(t) and u2(t) are the required 2 solutions 
           u1(a)=alpha1
           u2(a)=alpha2
           a<=t<=b
           a+h*n=b
    """
    tn=a
    y1i=alpha1
    y2i=alpha2
    y1s=[]
    y2s=[]
    y1s.append(alpha1)
    y2s.append(alpha2)
    y1i=y1s[0]
    y2i=y2s[0]
    for i in range(1,n+1):
        k1_1=f1(tn,y1i,y2i)
        k1_2=f2(tn,y1i,y2i)
        k2_1=f1(tn+(h*0.5),y1i+(h*k1_1*0.5),y2i+(h*k1_2*0.5))
        k2_2=f2(tn+(h*0.5),y1i+(h*k1_1*0.5),y2i+(h*k1_2*0.5))
        k3_1=f1(tn+(h*0.5),y1i+(h*k2_1*0.5),y2i+(h*k2_2*0.5))
        k3_2=f2(tn+(h*0.5),y1i+(h*k2_1*0.5),y2i+(h*k2_2*0.5))
        k4_1=f1(tn+h,y1i+h*k3_1,y2i+h*k3_2)
        k4_2=f2(tn+h,y1i+h*k3_1,y2i+h*k3_2)
        y1i1=y1i+(h*(1/6))*(k1_1+2*k2_1+2*k3_1+k4_1)
        y2i1=y2i+(h*(1/6))*(k1_2+2*k2_2+2*k3_2+k4_2)
        y1s.append(y1i1)
        y2s.append(y2i1)
        tn+=h
        y1i=y1i1
        y2i=y2i1

    return (y1s,y2s)        
        

def System2eq_RK2(alpha1,alpha2,f1,f2,h,n,a):
    """ 
    This functions solves an IVP system of 2 equations and returns the value at all the mesh
    points of both the equations, with Runge-Kutta 2.
    PROBLEM:
    Given: f1(t,u1,u2)=u1'
           f2(t,u1,u2)=u2'
           u1(t) and u2(t) is the required 2 solutions 
           u1(a)=alpha1
           u2(a)=alpha2
           a<=t<=b
           a+h*n=b
    """
    tn=a
    y1i=alpha1
    y2i=alpha2
    y1s=[]
    y2s=[]
    y1s.append(alpha1)
    y2s.append(alpha2)
    y1i=y1s[0]
    y2i=y2s[0]
    for i in range(1,n+1):
        k1_1=f1(tn,y1i,y2i)
        k1_2=f2(tn,y1i,y2i)
        k2_1=f1(tn+(h*0.5),y1i+(h*k1_1),y2i+(h*k1_2))
        k2_2=f2(tn+(h*0.5),y1i+(h*k1_1),y2i+(h*k1_2))
        y1i1=y1i+(h*(1/2))*(k1_1+k2_1)
        y2i1=y2i+(h*(1/2))*(k1_2+k2_2)
        y1s.append(y1i1)
        y2s.append(y2i1)
        tn+=h
        y1i=y1i1
        y2i=y2i1

    return (y1s,y2s)
