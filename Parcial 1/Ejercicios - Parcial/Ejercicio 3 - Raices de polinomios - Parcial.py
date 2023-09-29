#Calcular todas las raíces reales de f(x) = 3x^5 + 5x^4 -x^3
import numpy as np
import matplotlib.pyplot as plt

def Función (x=float):
    return 3*(x**5) + 5*(x**4) - (x**3)

x = np.linspace(-2,1,50)
y = Función(x)

plt.plot(x,y)
plt.axhline(y = 0,color='orange')
plt.show()

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

def Todas_las_Raíces (x, decimales=5):
    raíces = np.array([])
    for i in x: 
        una_raíz = Newton_Raphson(Función,Derivada,i)
        if una_raíz != False:
            Raíz = np.round(una_raíz, decimales)
            if Raíz not in raíces:
                raíces = np.append(raíces,Raíz)         
    raíces.sort()
    
    return raíces

raíces = Todas_las_Raíces(x)
print("las raices del polinomio son: " +str(raíces))