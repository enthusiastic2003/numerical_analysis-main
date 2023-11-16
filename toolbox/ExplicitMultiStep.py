import sys
sys.path.append('./toolbox/')

def AdamBashforthO2(df,x0,y0,y1,h,n):
    x1=x0+h
    y=0
    for i in range(2,n+1):
        y=y1+(h*0.5)*(3*(df(x1,y1))-df(x0,y0))
        x0=x1
        x1+=h
        y0=y1
        y1=y
    return y1


def AdamBashforthO3(df,x0,y0,y1,y2,h,n):
    x1=x0+h
    x2=x1+h
    y=0
    for i in range(3,n+1):
        y=y2+(h*(1/12))*(23*(df(x2,y2))-16*df(x1,y1)+5*df(x0,y0))
        y0=y1
        y1=y2
        y2=y
        x0=x1
        x1=x2
        x2=x2+h
    return y2

def AdamBashforthO4(df,x0,y0,y1,y2,y3,h,n):
    x1=x0+h
    x2=x1+h
    x3=x2+h
    y=0
    for i in range(4,n+1):
        y=y3+(h/24)*(55*df(x3,y3)-59*(x2,y2)+37*df(x1,y1)-9*(df(x0,y0)))
        y0=y1
        y1=y2
        y2=y3
        y3=y
        x0=x1
        x1=x2
        x2=x3
        x3=x3+h
        
    return y3
