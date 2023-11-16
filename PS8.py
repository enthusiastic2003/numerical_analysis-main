import math as m
import toolbox.ImplicitMultistep as IMS
import toolbox.RungeKutta as RK
import toolbox.ExplicitMultiStep as EMS
import toolbox.SystemSolver as SS

def df1(x,y):
    return ((x*m.exp(3*x))-2*y)

def f1(x):
    return (((1/5)*x*m.exp(3*x))-((1/25)*m.exp(3*x))+((1/25)*m.exp(-2*x)))

def df2(x,y):
    return 1+((x-y)**2)

def f2(x):
    return x+(1/(1-x))

def df3(x,y):
    return (1+(y/x))

def f3(x):
    return (x*m.log(x)+2*x)

def df4(x,y):
    return (y-x**2+1)

def p5f1(t,vu1,vu2):
    return (3*vu1)+(2*vu2)-((2*(t**2)+1)*m.exp(2*t))

def p5f2(t,vu1,vu2):
    return (4*vu1)+vu2+(((t**2)+(2*t)-4)*m.exp(2*t))

def u1(t):
    return ((1/3)*m.exp(5*t))-((1/3)*m.exp(-1*t))+m.exp(2*t)

def u2(t):
    return ((1/3)*m.exp(5*t))+((2/3)*m.exp(-1*t))+((t**2)*m.exp(2*t))

if __name__=="__main__":
    y1=RK.RK2_modified_euler(df1,0,0,0.2,1)
    yf=EMS.AdamBashforthO2(df1,0,0,y1,0.2,5)
    yfac=f1(1)
    print("Q1::a::actual=",yfac)
    print("ADAM_BASHFORTH_2=",yf) 
    print("DIFF=",abs(yfac-yf))
    print("#####################################")
    y1=RK.RK2_modified_euler(df2,2,1,0.2,1)
    yf=EMS.AdamBashforthO2(df2,2,1,y1,0.2,5)
    yfac=f2(3)
    print("Q1::b::actual=",yfac)
    print("ADAM_BASHFORTH_2=",yf) 
    print("DIFF=",abs(yfac-yf))
    print("#####################################")
    y1=RK.RK2_modified_euler(df3,1,2,0.2,1)
    yf=EMS.AdamBashforthO2(df3,1,2,y1,0.2,5)
    yfac=f3(2)
    print("Q1::c::actual=",yfac)
    print("ADAM_BASHFORTH_2=",yf) 
    print("DIFF=",abs(yfac-yf))
    print("#####################################")
##BEGINS MOULTON

    y1=RK.RK2_modified_euler(df1,0,0,0.2,1)
    yf=IMS.AdamsMoultonO2(df1,0,0,y1,0.2,5)
    yfac=f1(1)
    print("Q2::a::actual=",yfac)
    print("ADAM_MOULTON_2=",yf) 
    print("DIFF=",abs(yfac-yf))
    print("#####################################")

    y1=RK.RK2_modified_euler(df2,2,1,0.2,1)
    yf=IMS.AdamsMoultonO2(df2,2,1,y1,0.2,5)
    yfac=f2(3)
    print("Q2::b::actual=",yfac)
    print("ADAM_MOULTON_2=",yf) 
    print("DIFF=",abs(yfac-yf))
    print("#####################################")

    y1=RK.RK2_modified_euler(df3,1,2,0.2,1)
    yf=IMS.AdamsMoultonO2(df3,1,2,y1,0.2,5)
    yfac=f3(2)
    print("Q2::c::actual=",yfac)
    print("ADAM_MOULTON_2=",yf) 
    print("DIFF=",abs(yfac-yf))
    print("#####################################")

##BEGINS BASHFORTH 3
    y1=RK.RK2_modified_euler(df1,0,0,0.2,1)
    y2=RK.RK2_modified_euler(df1,0,0,0.2,2)
    yf=EMS.AdamBashforthO3(df1,0,0,y1,y2,0.2,5)
    yfac=f1(1)
    print("Q3::a1::actual=",yfac)
    print("ADAM_BASHFORTH_3=",yf) 
    print("DIFF=",abs(yfac-yf))
    print("#####################################")
    y1=RK.RK2_modified_euler(df2,2,1,0.2,1)
    y2=RK.RK2_modified_euler(df2,2,1,0.2,2)
    yf=EMS.AdamBashforthO3(df2,2,1,y1,y2,0.2,5)
    yfac=f2(3)
    print("Q3::b1::actual=",yfac)
    print("ADAM_BASHFORTH_3=",yf) 
    print("DIFF=",abs(yfac-yf))
    print("#####################################")
    y1=RK.RK2_modified_euler(df3,1,2,0.2,1)
    y2=RK.RK2_modified_euler(df3,1,2,0.2,2)
    yf=EMS.AdamBashforthO3(df3,1,2,y1,y2,0.2,5)
    yfac=f3(2)
    print("Q3::c1::actual=",yfac)
    print("ADAM_BASHFORTH_3=",yf) 
    print("DIFF=",abs(yfac-yf))
    print("#####################################")

############TODO:4############

######BEGINS P5###############
    (y1s,y2s)=SS.System2eq_RK4(1,1,p5f1,p5f2,0.2,5,0)
    y1ac=u1(1)
    y2ac=u2(1)
    print("Predicted u1(1)= ",y1s)
    print("predicted u2(1)=",y2s)
    print("Actual u1(1)=",y1ac)
    print("Actual u2(1)=",y2ac)


