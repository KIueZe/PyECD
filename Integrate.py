import numpy as np 
from scipy.integrate import quad
import math

def gauss(sigma, E, Ek, Rk):
    cons = float(2.4572717053474 * math.pow(10, -3)) # 10, 38
    # cons = 1/22.96 * np.power(np.pi, -1/2) * np.power(10, 40)
    sig_re = float(1/sigma)
    # reciprocal of sigma (half of the width of the CD band at 1/e peak height)
    expower = float(-math.pow((E-Ek) * sig_re, 2))

    epsilon = float(cons * sig_re * Ek * Rk * math.pow(math.e, expower))
    return epsilon

# integrate
val1,_=quad(lambda x:gauss(0.4, x, 6, 30) if gauss(0.4, x, 6, 30) > 0 else 0, # func
               5.5, # 积分下界
               6.5) # 积分上界
print ('积分结果：',_)
