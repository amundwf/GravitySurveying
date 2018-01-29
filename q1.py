import numpy as np
from matplotlib import pyplot as plt

#Interval I = [a,b]
a=0; b= 1
a0=1/3; b0=2/3
acc=100
d1=0.025; d2=0.25; d3=2.5

def F(x, a, b, d):
	return (b-x) / (d*(d**2+(x-b)**2)**(1/2)) - (a-x) / (d*(d**2+(x-a)**2)**(1/2))


def plotting():
	xvalues = np.linspace(a,b,acc)
	yvalues1 = np.zeros(acc)
	yvalues2 = np.zeros(acc)
	yvalues3 = np.zeros(acc)
	for i in range(acc):
		yvalues1[i] = np.log(F(xvalues[i], a0, b0, d1))
		yvalues2[i] = np.log(F(xvalues[i], a0, b0, d2))
		yvalues3[i] = np.log(F(xvalues[i], a0, b0, d3))
	plt.plot(xvalues, yvalues1, label='d1', lw=3)
	plt.plot(xvalues, yvalues2, label='d2', lw=3)
	plt.plot(xvalues, yvalues3, label='d3', lw=3)
	plt.legend()
	plt.xlabel("x", size=20)
	plt.ylabel("ln(F(x))", size=20)
	plt.grid()
	plt.show()

plotting()