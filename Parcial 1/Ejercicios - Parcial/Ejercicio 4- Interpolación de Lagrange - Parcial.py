import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sympy as sym
import math

def Lagrange(X, Y):

 x = sym.symbols('x')
 n = len(X)
 lagrange_polynomial = 0

 for i in range(n):
        term = Y[i]
        for j in range(n):
            if i != j:
                term *=(x - X[j])/(X[i] - X[j])
        lagrange_polynomial += term

 return sym.expand(lagrange_polynomial)

X = [1.4, 3.5, 5.6]
Y = [0.4007954931819738, 0.594128102489774, 0.29802795523938164]
print(Lagrange(X, Y))

def devcentral(x, h):
    min=10
    z=0
    i=0
    derivadas=[]
    d=((-0.0554912422401579*(x+h)**2 + 0.363970234266202*(x+h) + 5.55111512312578e-17)-(-0.0554912422401579*(x-h)**2 + 0.363970234266202*(x-h) + 5.55111512312578e-17))/(2*h)
    
    for i in d:
        derivadas.append(i)
    
   
    i=0  
    for derivada in derivadas:
        if round(derivada, 2)==0:
            indice=i
                
        i+=1
    
    distancia=x[indice]  
   
    return distancia
x = np.linspace(0,8,100)
distancia=devcentral(x, 1e-7)

ymax= -0.0554912422401579*(distancia)**2 + 0.363970234266202*(distancia) + 5.55111512312578e-17
Viy=(2*(9.8)*ymax)**(1/2)
tiempo=(Viy)/9.8
Vix=distancia/tiempo 
Vo=((Vix)**2+(Viy)**2)**(1/2)
print("El vector Inicial, que está definido por su magnitud y dirección es: ")
print("Magnitud: " +str(Vo) + " m/s")
Angulo=math.atan(Viy/Vix)
Angulo= Angulo*(180/math.pi)
print("Ángulo: " +str(Angulo) + " grados")


    


