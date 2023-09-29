from IPython.display import display, HTML
display(HTML("<style>.container { width:100% !important; }</style>"))

import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
from scipy import integrate
import math 
from mpl_toolkits.mplot3d import axes3d
sym.init_printing(use_unicode=True)

x = sym.Symbol('x',real=True)
y = sym.Symbol('y',real=True)

#Definición de varaibles para el cálculo de los límites de integración para el numerador
h=6.626e-34
Kb=1.3806e-23
c=3e+8
T=5772
lambda_1=4e-7
lambda_0=1e-7

V1=c/lambda_1
V0=c/lambda_0

#Integral del denominador utilizando cuadratura de Laguerre
n = 20
Roots,Weights = np.polynomial.laguerre.laggauss(n)

f= lambda x: (x**3/(math.e**(x)-1))/(math.e**(-x))
I_dem = 0
for i in range(n):
    I_dem += Weights[i]*f(Roots[i])



#Integral del numerador utilizando cuadratura de Legendre:
n = 20
a=(h*V1)/(Kb*T)
b=(h*V0)/(Kb*T)
def Integral(n, b, a):
    f1= lambda x: (x**3/(math.e**(x)-1))
    I_num = 0
    Roots, Weights = np.polynomial.legendre.leggauss(n)
    for i in range(n):
        I_num +=Weights[i]*f1(Roots[i]*(b-a)/2 + (b+a)/2)
    I_num = I_num*(b-a)/2
    
    return I_num

I_num=Integral(n, b, a)
f=(I_num/I_dem)*100
print("la fracción de rayos UV en términos de porcentaje es: "+str(f))



