import sys
sys.path.append('./toolbox/')

import RungeKutta as RK

def AdamsMoultonFixedPointO2(df,x0,y0,x1,y1,x2,y2_approx,h):
    y2=y2_approx
    for i in range(100):
        y=y1+(h/12)*(5*df(x2,y2)+8*df(x1,y1)-df(x0,y0))
        y2=y
    return y2
        
def AdamsMoultonFixedPointO3(df,x0,y0,x1,y1,x2,y2,x3,y3_pred,h):
    y3=y3_pred
    for i in range(100):
        y=y2+(h/24)*((9*df(x3,y3))+(19*df(x2,y2))-(5*(df(x1,y1)))+df(x0,y0))
        y3=y
    return y3

def AdamsMoultonO2(df,x0,y0,y1,h,n):
    x1=x0+h
    x2=x1+h
    y2_pred=RK.RK2_modified_euler(df,x1,y1,h/20,20)
    for i in range(n-1):
        y2=AdamsMoultonFixedPointO2(df,x0,y0,x1,y1,x2,y2_pred,h)
        x0=x1
        x1=x2
        x2+=h
        y0=y1
        y1=y2
        y2_pred=RK.RK2_modified_euler(df,x1,y1,h/20,20)
    return y2

def AdamsMoulton03(df,x0,y0,y1,y2,h,n):
    x1=x0+h
    x2=x1+h
    x3=x2+h
    y3_pred=RK.RK2_modified_euler(df,x2,y2,h/20,20)
    for i in range(n):
        y3=AdamsMoultonFixedPointO3(df,x0,y0,x1,y1,x2,y2,x3,y3_pred,h)
        x0=x1
        x1=x2
        x2=x3
        x3+=h
        y0=y1
        y1=y2
        y2=y3
        y3_pred=RK.RK2_modified_euler(df,x2,y2,h/20,20)
    return y3




    
