import SystemSolver as SS

def Deg2RK4(ypp,a,alpha,alphap,h,n):
    """ 
    This functions solves a Degree 2 IVP, and returns the function value and its derivative at all the mesh
    points , using RK4.
    PROBLEM:
    Given: ypp(t,y,yp) ;
           y(a)=alpha ;
           yp(a)=alphap ;
           
           a+h*n=b
    """
    yi=alpha 
    yip=alphap 
    yis=[yi] 
    yips=[yip]
    ti=a
    for i in range(1,n+1):
        k1_2=ypp(ti,yi,yip)
        k1_1=yip
        
        k2_2=ypp(ti+(0.5*h),yi+(h*k1_1*0.5),yip+(h*k1_2*0.5))
        k2_1=yip+(h*k1_2*0.5)

        k3_2=ypp(ti+(0.5*h),yi+(h*k2_1*0.5),yip+(h*k2_2*0.5))
        k3_1=yip+(h*k2_2*0.5)

        k4_1=ypp(ti+h,yi+h*k3_1,yip+h*k3_2)
        k4_2=yip+(h*k3_2)

        yi1=yi+(h/6)*(k1_1+(2*k2_1)+(2*k3_1)+k4_1)
        yip1=yip+(h/6)*(k1_2+(2*k2_2)+(2*k3_2)+k4_2)

        yis.append(yi1)
        yips.append(yip1)

        yi=yi1
        yip=yip1

        ti+=h

    return (yis,yips)