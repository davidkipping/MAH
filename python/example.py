import MAH
import simpleread
import math

# === SECTION 1: BASIC EXECUTION ===
# Here, I show a simple example of calling the MAH module

# = Read in the data ===
x = simpleread.simple_read("example_data.dat",2)
n = len(x)
MP=[]
RP=[]
for i in range(n):
  MP.append(x[i][0])
  RP.append(x[i][1])

# = Call the MAH subroutine =
result = MAH.MAH_(RP,MP)
RH2O = [[0 for k in range(2)] for i in range(n)]
RH2O = result[0]
RMAH = [[0 for k in range(2)] for i in range(n)]
RMAH = result[1]
PRMAH = [0 for k in range(2)]
PRMAH = result[2]
nvalid = [0 for k in range(2)]
nvalid = result[3]
valid=[[0 for k in range(2)] for i in range(n)]
valid = result[4]

# = Output the results==
zpure = [[0 for k in range(3)] for i in range(n)]
for i in range(n):
  zpure[i][0] = RH2O[i][0]
  zpure[i][1] = RMAH[i][0]
  zpure[i][2] = valid[i][0]
with open("100_H2O.dat", 'w') as f:
    f.writelines(' '.join(str(j) for j in i) + '\n' for i in zpure)

zdirty = [[0 for k in range(3)] for i in range(n)]
for i in range(n):
  zdirty[i][0] = RH2O[i][1]
  zdirty[i][1] = RMAH[i][1]
  zdirty[i][2] = valid[i][1]
with open("75_H2O_25_MgSiO3.dat", 'w') as f:
    f.writelines(' '.join(str(j) for j in i) + '\n' for i in zdirty)

# === SECTION 2: SOME SIMPLE STATISTICS ===
# Here, a few simple statistics are calculated using the results. These may
# be useful when preparing a paper.

# = Calculate median, negian and posian of RH2O & RMAH =
RH2O_sol = [[0 for j in range(3)] for k in range(2)]
RMAH_sol = [[0 for j in range(3)] for k in range(2)]
for k in range(2):
  result = MAH.best_values(RH2O,RMAH,nvalid,valid,k)
  for j in range(3):
    RH2O_sol[k][j] = result[0][j]
    RMAH_sol[k][j] = result[1][j]

# = Calculate median, negian and posian of RP & MP =
MP_sol = [0 for k in range(2)]
RP_sol = [0 for k in range(2)]
MP_sol = MAH.marginalize( MP )
RP_sol = MAH.marginalize( RP )

# = Summarize the results =
print '=== Observations ==='
print '$M_P =$ $',MP_sol[0],'_{-',MP_sol[1],'}^{+',MP_sol[2],'}$'
print '$R_P =$ $',RP_sol[0],'_{-',RP_sol[1],'}^{+',RP_sol[2],'}$'
print ' '
for k in range(2):
  if k == 0:
    print '=== 100%-H2O model ==='
  else:
    print '=== 75%-H2O-25%-MgSiO3 model ==='
  print '$R_{H2O} =$ $',RH2O_sol[k][0],'_{-',RH2O_sol[k][1],'}^{+',\
               RH2O_sol[k][2],'}$'
  print '$R_{MAH} =$ $',RMAH_sol[k][0],'_{-',RMAH_sol[k][1],'}^{+',\
               RMAH_sol[k][2],'}$'
  print '$P(R_{MAH}>0) =$ ',PRMAH[k]*1e2,'\% (',\
               math.sqrt(2.0)*MAH.inverf(PRMAH[k]),'-sigma)'
  print ' ' 
