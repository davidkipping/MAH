import math

# ==============================================================================
def MAH_(RP,MP):
  n = len(RP)
  # Compute RMAH
  nvalid=[0 for k in range(2)]
  valid=[[0 for k in range(2)] for i in range(n)]
  RH2O=[[0 for k in range(2)] for i in range(n)]
  RMAH=[[0 for k in range(2)] for i in range(n)]
  temp=[0 for k in range(2)]
  for k in range(2):
    for i in range(n):
      temp = R_water(MP[i],k)
      RH2O[i][k] = temp[0]
      valid[i][k] = temp[1]
      if valid[i][k] == 1:
        nvalid[k] = nvalid[k] + 1 # Increase the count of valid points
      RMAH[i][k] = RP[i] - RH2O[i][k]
  # What percentage of inputs were valid?
  fvalid = []
  for k in range(2):
    fvalid.append( float(nvalid[k])/float(n) )
    #print 'For model ',k,' ',fvalid[k]*1e2,'% of trials valid'
  # Compute P(RMAH>0)
  npositive=[0 for k in range(2)]
  PRMAH=[0 for k in range(2)]
  for k in range(2):
    for i in range(n):
      if valid[i][k] == 1 and RMAH[i][k] > 0:
        npositive[k] = npositive[k] + 1
    PRMAH[k] = float(npositive[k])/float(nvalid[k])
    #print 'For model ',k,' P(RMAH>0) = ',PRMAH[k]
  # Return the result
  result = [ RH2O ]
  result.append( RMAH )
  result.append( PRMAH )
  result.append( nvalid )
  result.append( valid )
  return result
# ==============================================================================

# ==============================================================================
def R_water(M,flag):
  # Define various constants
  p = [+1.409482429344076e0,\
       +3.942477439199774e-1,\
       +5.014920562075038e-2,\
       +2.513336445152832e-3,\
       -4.557149516701664e-4,\
       -9.717343773829619e-5,\
       -3.900230062481005e-6,\
       +1.776907877676882e-7]
  q = [+1.346431690428415e+0,\
       +3.797476623436896e-1,\
       +4.669407861902054e-2,\
       +1.992173679766138e-3,\
       -3.468932022678234e-4,\
       -7.637885099727161e-5,\
       -6.314645065750879e-6,\
       -1.981417530222748e-7]
  MPmin = [4.854443348229839e-4,6.196844780385370e-5]
  MPmax = [4.864096638526573e+2,3.393232615628230e+2]
  a=[[0 for k in range(2)] for i in range(len(p))]
  for j in range(len(p)):
    a[j][0] = p[j]
    a[j][1] = q[j]
  # Check if valid
  if M > MPmin[flag] and M < MPmax[flag]:
    valid = 1 # M is inside range of Zeng & Sasselov (2013) models
  else:
    valid = 0 # M is outside range of Zeng & Sasselov (2013) models
  # Calculate Rwater
  logM = math.log(M)
  Rwater = a[0][flag]
  for j in range(7):
    Rwater = Rwater + a[j+1][flag]*(logM)**(j+1)
  # Return result
  result = [ Rwater, valid ]
  return result
# ==============================================================================

# ==============================================================================
def sort_array(x):
  y = sorted(x)
  return y
# ==============================================================================

# ==============================================================================
def best_values(RH2O,RMAH,nvalid,valid,k):
  # Define a list of valid results only
  n = len(RH2O)
  RH2O_temp = [0 for m in range(nvalid[k])]
  RMAH_temp = [0 for m in range(nvalid[k])]
  m = 0
  for i in range(n):
    if valid[i][k] == 1:
      m = m + 1
      RH2O_temp[m-1] = RH2O[i][k]
      RMAH_temp[m-1] = RMAH[i][k]
  RH20_sol = marginalize(RH2O_temp)
  RMAH_sol = marginalize(RMAH_temp)
  result = [ RH20_sol ]
  result.append( RMAH_sol )
  return result
# ==============================================================================

# ==============================================================================
def marginalize(x):
  # Useful constants
  n = len(x)
  negian_I = int( math.floor( 0.1586552539314571*n ) )
  median_I = int( round( 0.5*n ) )
  posian_I = int( math.ceil( 0.8413447460685429*n ) )
  # Sort the list
  x_sort = sort_array(x)
  # Extract negian, median and posian
  x_sol = [0 for i in range(3)]
  x_sol[0] = x_sort[median_I]
  x_sol[1] = x_sort[median_I] - x_sort[negian_I]
  x_sol[2] = x_sort[posian_I] - x_sort[median_I]
  return x_sol
# ==============================================================================

# ==============================================================================
def inverf(x):
  # Constants
  a_winitzki = 0.1400122886866666
  b_winitzki = 4.546884979448284
  # Compute factor
  if x<0.0:
    factor = -1.0
  else:
    factor = 1.0
  # Compute LDCs
  logxsq = math.log(1.0-x*x)
  y = b_winitzki + 0.5*logxsq
  y = math.sqrt( y*y - (logxsq/a_winitzki) ) - y
  y = factor*math.sqrt(y)
  return y
# ==============================================================================
