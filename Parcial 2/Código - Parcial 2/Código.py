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

"Según la página del IDEAM, el porcentaje de ultravioleta es de 7.2%, mientras que la calculada es de 12% aproximadamente."
"Esto se debe a que la radiación solar emitida por el sol (12% aprox) al pasar por la atmósfera, sufre un proceso de"
"debilitamiento por la difusión, reflexión en las nubes y de absorción por las moléculas de gases"
"(como el ozono y el vapor de agua) y por partículas en suspensión, la radiación solar alcanza la superficie"
"terrestre oceánica y continental que la refleja o la absorbe"


#Parte teórica. 
#Punto A y B
x = sym.Symbol('x',real=True)
y = sym.Symbol('y',real=True)


def GetLaguerreRecursive(n,x):
    if n==0:
        poly = sym.Number(1)
    elif n==1:
        poly = 1 - x
    else:
        poly = ((2*(n-1)+1-x)*GetLaguerreRecursive(n-1,x)-(n-1)*GetLaguerreRecursive(n-2,x))/n
   
    return sym.expand(poly,x)

def GetDLaguerre(n,x):
    Pn = GetLaguerreRecursive(n,x)
    return sym.diff(Pn,x,1)

def GetNewton(f,df,xn,itmax=10000,precision=1e-5):
    
    error = 1.
    it = 0
    
    while error >= precision and it < itmax:
        
        try:
            
            xn1 = xn - f(xn)/df(xn)
            
            error = np.abs(f(xn)/df(xn))
            
        except ZeroDivisionError:
            print('Zero Division')
            
        xn = xn1
        it += 1
        
    if it == itmax:
        return False
    else:
        return xn
    
def GetRoots(f,df,x,tolerancia = 5):
    
    Roots = np.array([])
    
    for i in x:
        
        root = GetNewton(f,df,i)

        if  type(root)!=bool:
            croot = np.round( root, tolerancia )
            
            if croot not in Roots:
                Roots = np.append(Roots, croot)
                
    Roots.sort()
    
    return Roots


def GetAllRootsGLag(n):
    xn = np.linspace(0,n + (n-1)*math.sqrt(n), 1000)
    
    Laguerre = []
    DLaguerre = []
    
    for i in range(n+1):
        Laguerre.append(GetLaguerreRecursive(i,x))
        DLaguerre.append(GetDLaguerre(i,x))
    
    poly = sym.lambdify([x],Laguerre[n],'numpy')
    Dpoly = sym.lambdify([x],DLaguerre[n],'numpy')
    Roots = GetRoots(poly,Dpoly,xn)

    if len(Roots) != n:
        ValueError('El número de raíces debe ser igual al n del polinomio.')
    
    return Roots

print("El polinomio de orden 2 es: ")
print(GetLaguerreRecursive(2, x))
print("Las dos raíces del segundo polinomio son")
print(GetAllRootsGLag(2))
x_l=GetAllRootsGLag(2)


x = sym.symbols('x')
f1= lambda x: (x-3.41421)/(0.58579 -3.41421)
f2=lambda x:(-x)*(x-0.58579)/(3.41421-0.58579)

def Laguerre_1(n):
    Roots,Weights=np.polynomial.laguerre.laggauss(n)
    suma= np.sum(Weights*f1(Roots))
    return suma

primer_peso=Laguerre_1(10)
print(primer_peso)

def Laguerre_2(n):
    Roots,Weights=np.polynomial.laguerre.laggauss(n)
    suma= np.sum(Weights*f2(Roots))
    return suma

segundo_peso=Laguerre_2(10)
print(segundo_peso)

f3=lambda x: x**3

def Laguerre_3(n):
    Roots,Weights=np.polynomial.laguerre.laggauss(n)
    suma= np.sum(Weights*f3(Roots))
    return suma

ultimo=Laguerre_3(20)
print(ultimo)










