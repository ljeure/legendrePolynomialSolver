e1 = 1e-6
e2 = 1e-6

def P(n, x):
  
  if n == 0:
    return 1
  if n == 1:
    return x
  else:
    c = 2.0*(n-2) + 3.0
    b = 1.0*(n-2) + 2.0
    a = 1.0*(n-2) + 1.0
    value = c/b*x*P(n-1, x) - a/b*P(n-2,x)
    return value

# not currently used
def F(m, x):
  if x == 0:
    return 0
  if m == 1:
    return (x**2 - 1)/(2*x)

  return 1/(2*m*(x - (m - 1)*F(m-1,x)))

# P'/P
def S1(m, x):
  num = m*x - m*P(m-1,x)/P(m,x)
  denom = x**2 - 1
  return num/denom

def S2(m, x):
  num = m*(m+1) + S1(m,x) * ((1-x**2)* S1(m,x) - 2*x)
  denom = 1-x**2
  return -num/denom

def findLegendreRoots(n):

  roots=list()
  converged=list()
  s1_tilde=list()
  s2_tilde=list()

  # set guesses with log scale
  for i in range(n/2):
    roots.append(- 2**(-.5*(i+1)) +1)
    converged.append(False)
    s1_tilde.append(0)
    s2_tilde.append(0)

  for i in range(n/2):
    roots.append(-roots[i])
    converged.append(False)
    s1_tilde.append(0)
    s2_tilde.append(0)

  if n % 2 == 1:
    roots.append(0.0)
    converged.append(True)
    s1_tilde.append(0)
    s2_tilde.append(0)

  all_roots_converged = False

  roun = 0

  # use the Alberth-Housholder_n method to nudge guesses towards roots
  while not all_roots_converged:
    roun += 1
    
    # set S tildes
    for i in range(n):
      if not converged[i]:
        sum1 = 0
        sum2 = 0
        for j in range(n):
          if j != i:
            sum1 += 1/(roots[i] - roots[j])
            sum2 += -1/(roots[i] - roots[j])**2
    
        s1_tilde[i] = S1(n, roots[i]) - sum1
        s2_tilde[i] = S2(n, roots[i]) - sum2

        # nudge roots towards actual roots
        #for i in range(n):
        #if not converged[i]:

        # householder method 1: Newton-Raphson
        #u_new = roots[i] - 1/(s1_tilde[i])
        
        # householder method 2  Halley    
        u_new = roots[i] - 2*s1_tilde[i] / (s1_tilde[i]**2 - s2_tilde[i])

        # if this is the actual root
        if abs(u_new - roots[i]) < e1:
          if P(n, u_new) < e2:
            converged[i] = True
        roots[i] = u_new

    for i in range(n):
      all_roots_converged = converged[i]
      if not all_roots_converged:
        break

  return sorted(roots)


# calculate the weights to be used for Gauss Legendre quadrature
def calculateWeights(roots):
  n = len(roots)
  weights = list()
  for i in range(n):
    weights.append( - (2*roots[i]**2 - 2) / (n*P(n-1, roots[i]))**2)
  
  return weights

lRoots = findLegendreRoots(15)
print 
print "roots:"
print lRoots
print
print "weights:"
print calculateWeights(lRoots)

