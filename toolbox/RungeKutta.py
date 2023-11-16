import sys
sys.path.append('./toolbox/')

def RK2_modified_euler(f,x0,y0,h,n):
  x=x0
  y=y0
  for i in range(n):
    k1=h*f(x,y)
    k2=h*f(x+h,y+k1)
    y+=(k1+k2)/2
    x+=h
  return y

def RK4_modified_euler(f,x0,y0,h,n):
  x=x0
  y=y0
  for i in range(n):
    k1=h*f(x,y)
    k2=h*f(x+(h/2),y+(1/2)*k1)
    k3=h*f(x+(h/2),y+(1/2)*k2)
    k4=h*f(x+(h),y+k3)
    y=y+(1/6)*(k1+2*k2+2*k3+k4)
    x=x+h
  return y


