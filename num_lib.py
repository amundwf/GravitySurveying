def bisection(a,b,function, tol, const):   #Returns x: functionx(x) = const. [a,b] is the interval, tol is the accuracy
	while (b-a)/2 > tol:
		c = (a+b)/2
		if function(c) == 0:
			return c
		if (function(c)-const)*(function(a)-const) < 0:
			b = c
		else:
			a = c
	return (a+b)/2

#Fixed point interation algoritm, returns 0,0 if iterations >= steps
def fpi(g, x0, tol, steps): #returns x,S: g(x)=x to the order of tol for initial x0
	x = x0
	counter = 0
	while True:
		counter += 1
		y = g(x)
		z = g(y)
		error1 = abs(x-y)
		error2 = abs(y-z)
		error = (error1 + error2)/2
		S = error2/error1
		if error < tol:
			break
		if counter >= steps:
			return 0, 0
		x = y
	return x, S
