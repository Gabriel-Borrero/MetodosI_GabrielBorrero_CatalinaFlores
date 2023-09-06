import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sympy as sym
import math

def Función (x=float):
    return math.e**(-x) -x

x = np.linspace(-2,1,50)
y = Función(x)

plt.plot(x,y)
plt.axhline(y = 0,color='orange')
#plt.show()

def Derivada(f,x,h=1e-6):
    return (f(x+h)-f(x-h))/(2*h)

def Newton_Raphson(f,derivada_f,xn,freno_iteraciones=100,precision=1e-8):
    error = 1.
    it = 0
    while error > precision and it < freno_iteraciones:
        try:  
            xn1 = xn - f(xn)/derivada_f(f,xn)
            # Criterio de parada
            error = np.abs(f(xn)/derivada_f(f,xn)) 
        except ZeroDivisionError:
            print('Division por cero')
        xn = xn1
        it += 1
    if it == freno_iteraciones:
        return False
    else:
        return xn
raíz = round(Newton_Raphson(Función, Derivada,1.), 5)

print("Para el punto a, La raíz es: "+str(raíz))

#Punto B:
x0=0
x1=2
y_positivo=math.e**(-x0) -(x0)
y_negativo=math.e**(-x1) -(x1)

#Comprobando el teorema de bolzano, se tiene que:

Bolzano=round(y_positivo*y_negativo,5)
print("Al multiplicar f(xo) y f(x1), obtenermos un valor que cumple con la condición de Bolzano de ser menor a cero, el cual es: "+str(Bolzano))

x2=(x0+x1)/2
print("Para el punto C, se tiene que x2 es igual a: "+str(x2))

#Para el punto D, primero calculamos el polinomio mediante interpolación de newton
def InterpolacionNewton(X,Y,x):
    
    sum_ = Y[0]
    
    Diff = np.zeros(( X.shape[0],Y.shape[0] ))
    h = X[1]-X[0]
    
    Diff[:,0] = Y

    poly = 1.
    
    for i in range(1,len(X)):
        
        poly *= (x-X[i-1])
        
        for j in range(i,len(X)):
            
            Diff[j,i] = Diff[j,i-1] - Diff[j-1,i-1] 
    
        sum_ += poly*Diff[i,i]/(math.factorial(i)*h**(i))
        
    return sym.expand(sum_)
    
x = sym.symbols('x')
X = np.array([0, 1, 2])
Y = np.array([1, -0.632, -1.865])

Polinomio=InterpolacionNewton(X, Y, x)
print("El polinomio calculado mediante interpolación de Newton es: " +str(Polinomio))

#Ahora, usando los coeficientes del polinomio antes calculado, encontraremos x3


a= 0.1995
b= -1.8315
c= 1.0
x3_suma= (-2*c)/(b+math.sqrt((b**2)-4*a*c))
x3_resta=(-2*c)/(b-math.sqrt((b**2)-4*a*c))

def coeficientes(x, X, Y, freno_iteraciones=100, precision=1e-10):
    a= 0.1995
    b= -1.8315
    c= 1.0
    x3_suma= (-2*c)/(b-math.sqrt((b**2)-4*a*c))
    fx3=math.e**(-x3_suma) -x3_suma
    X_nuevo=[]
    it = 0
    Y_nuevo=[]
    while  it < freno_iteraciones and abs(fx3)<precision:
        
        sum_ = Y[0]
        coeficientes=[]
        Diff = np.zeros(( X.shape[0],Y.shape[0] ))
        h = X[1]-X[0]
    
        Diff[:,0] = Y

        poly = 1.
    
        for i in range(1,len(X+1)):
        
            poly *= (x-X[i-1])
        
            for j in range(i,len(X+1)):
            
                Diff[j,i] = Diff[j,i-1] - Diff[j-1,i-1] 
        for i in range(0,len(X+1)):
            coeficientes.append(Diff[i,i])
        
            for i in coeficientes:
                a=coeficientes[2]
                b=coeficientes[1]
                c=coeficientes[0]
            
                if i==1:
                        if coeficientes[i]<0:
                             x3=(-2*c)/(b-math.sqrt((b**2)-4*a*c))
                        else:
                             x3=(-2*c)/(b+math.sqrt((b**2)-4*a*c))
        
        fx3= math.e**(x3) -x3
        it+=1
        for i in len(X):
            X_nuevo[i]=X[i+1]
            X.add(-1, x3)  
        for i in len(X):
            Y_nuevo[i]=Y[i+1]
            Y.add(-1, fx3)
            
        X=X_nuevo
        Y=Y_nuevo
        
              
    nuevo_polinomio=InterpolacionNewton(X,Y,x)
        
    return nuevo_polinomio


x = sym.symbols('x')
X = np.array([0, 1, 2])
Y = np.array([1, -0.632, -1.865])    
print(coeficientes(x, X, Y, freno_iteraciones=100, precision=1e-10)) 

nueva_raiz=Newton_Raphson(coeficientes(x, X, Y, freno_iteraciones=100, precision=1e-10),derivada_f,xn,freno_iteraciones=100,precision=1e-8)

print(nueva_raiz)
           
