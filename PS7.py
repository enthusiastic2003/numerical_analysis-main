import toolbox.RungeKutta as RK
import math as m

def df(x,y):
    return 1+((x-y)**2)
def f(x):
    return x+(1/(1-x))

print("AC=",f(3))
print("PD=",RK.RK2_modified_euler(df,2,1,0.1,10))