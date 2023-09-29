
import numpy as np
import matplotlib.pyplot as plt


x= np.arange(0.1,1.1,0.01)

tangente= lambda x: (np.tan(x))**0.5

dertan= lambda x: 1/(2*(np.tan(x))**0.5*np.cos(x)**2)

y=dertan(x)

def discretizacion(x,f,h=0.01):
    return (-3*f(x) + 4*f(x+h) - f(x+2*h))/(2*h)
progresiva= discretizacion(x,tangente)


def DerivativeC(x,f,h):
    return (f(x+h)-f(x-h))/(2*h)
central=DerivativeC(x,tangente,0.01)


fig= plt.figure(figsize=(10,10))
ax1=fig.add_subplot(2,2,1)
ax2=fig.add_subplot(2,2,2)
ax3=fig.add_subplot(2,2,3)
ax4=fig.add_subplot(2,2,4)


ax1.set_title("Derivada Central")
ax1.plot(x,central, "r*")
ax1.plot(x,y, "k")


ax2.set_title("Derivada Progresiva")
ax2.plot(x,progresiva, "g*")
ax2.plot(x,y, "k")


ax3.set_title("Error Derivada Central")
ax3.scatter(x,np.abs(y-central), color="r",s=0.5)


ax4.set_title("Error Derivada Progresiva")
ax4.scatter(x,np.abs(y-progresiva), color="g", s=0.5)
plt.show()

print("El valor máximo para el error en la derivada progresiva es: " + str(np.max(np.abs(y-progresiva))))
print("El valor máximo para el error en la derivada central es: " + str(np.max(np.abs(y-central))))


#CONCLUSIONES:

#En resumen, podemos concluir que tanto la derivada progresiva como la derivada central ofrecen 
#una aproximación muy precisa a la derivada real, hasta el punto en que, al representarlas gráficamente 
#junto a la derivada real, son prácticamente indistinguibles. En cuanto al error en los puntos nodales, 
#notamos que la derivada central es ligeramente más precisa que la derivada progresiva, 
#ya que el error máximo de la derivada progresiva es de 0.003, mientras que el de la derivada central es de 0.002. 
#No obstante, consideramos que la derivada progresiva calculada mediante el polinomio de interpolación 
#es bastante similar a la derivada central y a la derivada real. Por lo tanto, concluimos que es aceptable utilizar
#este polinomio de interpolación para estimar la derivada de una función.





