import numpy as np
import pickle
from TMA4320_invit.Project1_GravitySurveying.test_example import analytical_solution

a=0; b=1; omega = 3*np.pi; gamma=-2
d=0.025
Nc = Ns = 40
Nq = 40

def K(x, y, d):
	return d * (d**2+(y-x)**2)**(-3/2)

def fredholm_rhs(xc, F):
	# Expects collocation points xc, F(x) = function for the force in a coll.point
	# Returns the vector F = [F(x_0^c),  ..., F((x_(Ns-1)^c)]
	Nc = xc.shape[0]
	b = np.zeros(Nc)
	for i in range(len(b)):
		b[i] = F[(i,0)]
	return b
	
def fredholm_lhs(xc, xs, xq, w, K, Nq, d):
	# Expects collocation points xc, source points xs, quadrature points xq, quadrature weights w (All np.arrays)
		# and K = integral kernel
	# Returns The np matrix A s.t. A*rho = F
	Nc = xc.shape[0]; Ns = xs.shape[0]
	A = np.zeros((Nc, Ns))
	a = 0
	for i in range(Nc):
		for j in range(Ns):
			a = 0
			for k in range(Nq):
				a += w[k] * K(xc[k], xq[k], d) * lagrangePolyJ(j, xq[k], xs)
			#A[(i,j)] = a
			print(a)
	return A

def lagrangePolyJ(j, x, xs):
	# Expects the basis poly nr. j, the evaluation point x, the list of interpol. points xs of length >= j
	# Returns the value of the j-th Lagrange polynomial in x
	a = 1; b = 1
	for m in range(len(xs)):
		if m != j:
			a *= (x-xs[m])
			b*= (xs[j]-xs[m])
	return a/b

def rho(x):
	return np.sin(omega*x)*np.e**(gamma*x)

def conformalMapping(x, a, b):
	# Expects a point x E [-1,1] and interval limits a, b
	# Returns w in the mapping from x E [-1,1] to w E [a,b]
	return round(x*(b-a)/2 + (a+b)/2, 10)

def chebyshev(f, a, b, n):
	# Expects an integer n = the number of interpolation points (on the interval [-1,1], function f(x) to be
		# interpolated and interval limits a, b
	# Returns the np array [(x0,f(x0)), (x1,f(x1), ..., (x(n-1)f(x(n-1)))] = The interpolation points for the n-th
		# chebyshev polynomial
	A = np.zeros((n,2))
	for i in range(n):
		A[(i,0)] = conformalMapping(np.cos((np.pi/2 + np.pi*i)/n), a, b)
	for i in range(n):
		A[(i,1)] = f(A[(i,0)])
	return A

def midpointNC_x_w(f, a, b, N):
	# Expects a function f(x), interval [a,b] and N number of panels
	# Returns the np array [x0, x1, ..., x(n-1)] and the np array [w0, w1, ..., w(n-1)] for the x-values and weights
		# for a midpoint quadrature of f(x) on N panels.
	x = np.zeros(N); w = np.zeros(N)
	h = (b - a) / N
	for i in range(N):
		x[i] = a + h/2 + i*h
		w[i] = f(x[i])
	return x, w

xc = chebyshev(rho, a, b, Nc); xs = chebyshev(rho, a, b, Ns)
xq, w = midpointNC_x_w(rho, 0, 1, Nq)

f = pickle.load( open( "F.pkl", "rb" ))
F = f(xc, d)

b = fredholm_rhs(xc, F)
print(b)

A = fredholm_lhs(xc, xs, xq, w, K, Nq, d)
print(A)