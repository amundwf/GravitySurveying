import numpy as np
from matplotlib import pyplot as plt
import TMA4320_invit.Project1_GravitySurveying.proj1lib as lib

#Interval I = [a,b]
a=0
b=1
acc=100
d1=0.025; d2=0.25; d3=2.5

def rho(x):
	if (1/3 < x and x < 2/3):
		return 1
	else:
		return 0

def f(x,y,d):
	return d / ((d**2+(y-x)**2)**(3/2)) * rho(y)
def g(y):
	return

xvalues = np.linspace(a,b,acc)
yvalues = np.zeros(acc)

for i in range(acc):
	
	yvalues[i] = lib.simpson(f(xvalues[i], xvalues[i], d1), a, b, acc)
	print(xvalues[i], yvalues[i])

plt.plot(xvalues, yvalues)
plt.show()
