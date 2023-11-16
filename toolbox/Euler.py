import sys
sys.path.append('./toolbox/')

import RungeKutta as RK
def taylor_ivf(func,x0,y0,h,n):
  x=x0
  y=y0
  for i in range(n):
    y+=h*func(x,y)
    x+=h
  return y

def taylor_ivf_o2(df,ddf,x0,y0,h,n):
  x=x0
  y=y0

  for i in range(n):
    y+=h*df(x,y)+(h**2/2)*df(x,y)
    x+=h
  return y

def trap_iter(df,y0,x0,y1_approx,h):
  x1=x0+h
  y1=y1_approx
  for i in range(20):
    y1=y0+0.5*h*(df(x0,y0)+df(x1,y1))
  return y1

def trapezoid_method(df,x0,y0,h,n):
  x=x0
  y=y0
  w=RK.RK2_modified_euler(df,x,y,h,1)
  for i in range(n):
    y=trap_iter(df,y,x,w,h)
    x+=h
    w=RK.RK2_modified_euler(df,x,y,h,n)
  return y

