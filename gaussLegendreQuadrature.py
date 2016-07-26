import numpy as np

'''
a polynomial
'''
def F(x):
  return 3*x^4+2*x*2+5.3*x-2.7


'''
integrates a function F on the interval (a,b) using
Gause Legendre quadrature. n is the number of divisions
'''
def integrateF(a, b, n):
  n-=1
  if n>5:
    print "I need more quadrature to do that"
    exit()

  abscissa = [[0],[-0.5773502691896257,0.5773502691896257], \
      [-0.7745966692414834, 0, 0.7745966692414834], \
      [-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, \
      0.8611363115940526],
      [-0.9061798459386640, -0.5384693101056831, 0, 0.5384693101056831, \
          0.9061798459386640]]

  weight = [[2],[1, 1], [5./9., 8./9., 5./9.], \
      [0.3478548451374538, 0.6521451548625461, 0.6521451548625461, \
      0.3478548451374538], \
      [0.2369268850561891, 0.4786286704993665, 0.5688888888888889, \
      0.4786286704993665, 0.2369268850561891]]

  integral = 0

  for i in range(n+1):
    integral += (b-a)/2.* weight[n][i] * F((b-a)/2.*abscissa[n][i] +(b+a)/2.)

  return integral

print integrateF(-95,189,5)

