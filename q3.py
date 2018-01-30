import numpy as np

def fredholm_rhs(xc, F):
	#Expects collocation points xc, F(x) = function for the force in a coll.point
	#Returns the vector F = [F(x_0^c),  ..., F((x_(Ns-1)^c)]
	Nc = xc.shape[0]
	b = np.zeros(Nc)
	for i in range(len(b)):
		b[i] = F(xc[i])
	return b

def F(x):
	#Example function
	return 7

def fredholm_lhs(xc, xs, xq, w, K):
	#Expects collocation points xc, source points xs, quadrature points xq, quadrature weights w (All np.arrays)
		#and K = integral kernel
	#Returns The np matrix A s.t. A*rho = F
	Nc = xc.shape[0]; Ns = xs.shape[0]
	A = np.zeros((Nc, Ns))
	a = 0
	for i in range(Nc):
		for j in range(Ns):
			a = 0
			for k in range(nq):
				a += w[k] * K(xc[i], xq[k]) * lagrangePolyj(j, xq[k], xs)
			A[(i,j)] = a


def lagrangePolyJ(j, x, xs):
	#Expects the basis poly nr. j, the evaluation point x, the list of interpol. points xs of length >= j
	#Returns the value of the j-th Lagrange polynomial in x
	a = 1; b = 1
	for m in range(len(xs)):
		if m != j:
			a *= (x-xs[m])
			b*= (xs[j]-xs[m])
	return a/b