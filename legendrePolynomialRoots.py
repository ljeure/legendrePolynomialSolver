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

  if n%2 == 1:
    roots.append(0.)
    converged.append(True)
    s1_tilde.append(0)
    s2_tilde.append(0)


  all_roots_converged = False

  # use the Alberth-Housholder_n method to nudge guesses towards roots
  while not all_roots_converged:
    
    # set S tildes
    for i in range((n+1)/2):
      if not converged[i]:
        sum1 = 0
        sum2 = 0
        for j in range((n+1)/2):
          if j != i:
            sum1 += 1/(roots[i] - roots[j])
            sum2 += -1/(roots[i] - roots[j])**2
    
        s1_tilde[i] = S1(n, roots[i]) - sum1
        s2_tilde[i] = S2(n, roots[i]) - sum2

        # householder method 1: Newton-Raphson
        #u_new = roots[i] - 1/(s1_tilde[i])
        
        # householder method 2  Halley    
        u_new = roots[i] - 2*s1_tilde[i] / (s1_tilde[i]**2 - s2_tilde[i])
        u_old = roots[i]
        roots[i] = u_new

        # if this is the actual root
        if abs(u_new - u_old) < e1:
          if P(n, u_new) < e2:
            converged[i] = True

            # if this root equals another root
            for j in range((n+1)/2):
              if j != i:
                if abs(roots[j] - roots[i]) < e1:
                  print "root is not actually converged: "
                  print roots[i]
                  print

                  # reset the root to its original guess
                  roots[i] = - 2**(-.5*(i+1)) +1
                  converged[i] = False


    for i in range((n+1)/2):
      all_roots_converged = converged[i]
      if not all_roots_converged:
        break
  
  # add negative roots
  for i in range(n/2):
    roots.append(-roots[i])

  return roots

# calculate the weights to be used for Gauss Legendre quadrature
def calculateWeights(roots):
  n = len(roots)
  weights = list()
  for i in range(n):
    weights.append( - (2*roots[i]**2 - 2) / (n*P(n-1, roots[i]))**2)
  
  return weights


n = 64
lRoots = sorted(findLegendreRoots(n))
print
print "for n:"
print n
print 
print "roots:"
print lRoots
print
print "weights:"
print calculateWeights(lRoots)
print
