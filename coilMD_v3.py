import numpy as np
from openmmlib import Simulation
import joblib
import simtk.unit as units
import sys

filename = str(sys.argv[1])
block = joblib.load(filename)
#data = block
data = np.array(block['data'])
N = data.shape[0]

lineMonomersfile = str(sys.argv[2])
with open(lineMonomersfile) as f: lines = f.read().splitlines()
lmono = []
for i in range(len(lines)): lmono.append(int(lines[i]))
M = len(lmono)
dmono = data[lmono,:]
print(M)

# include condensins bindings
condensinsfile = str(sys.argv[3])
condensinIIdata = joblib.load(condensinsfile)
nLoopsII = condensinIIdata.shape[0]

outfolder = str(sys.argv[4]) 

rh = float(sys.argv[5])
pitch = 2*rh # in the case of the chromonema behaving as a rigid rope
print('The helical path radius is %d' % rh)
turnPath = np.sqrt((2*np.pi*rh)**2 + (2*rh)**2)
print('The turn path length is %d' % turnPath)

# the new positions of the chosen monomers depend on the position of the first monomer to be fixed
#firstMono = [0,0,0] # to put it in the origin, or any other coordinate
#firstMono = dmono[0,:] # to start with the first monomer in the list
# to start from the middle
firstMono = M//2
print(M//2,lmono[M//2],'first monomer to be fixed')
lmonof = lmono[M//2:]
lmonor = np.flip(lmono[:M//2+1])
print(len(lmonof),len(lmonor))

# r defines the direction of the straight line that M monomers will be pulled to
#r = np.array([1,0,0])
#r = (dmono[1,:]-dmono[0,:])/(np.sqrt(np.sum((dmono[1,:]-dmono[0,:])**2))) # if it starts from the beggining
# two directions if it starts from the middle
#rmod = np.sqrt(np.sum((dmono[firstMono+1,:]-dmono[firstMono-1,:])**2))
#rf = (dmono[firstMono+1,:]-dmono[firstMono-1,:])/rmod
#rr = -rf
#print (rf,rr)

# s is a list of the original distances between each chosen monomers
s = np.sqrt(np.sum(np.diff(dmono,axis=0)**2,1))
sums = np.sum(s)
print("The axis is %2.f nm long" % sums)
nTurns = sums/turnPath
print("The helix will have %2.f turns" % nTurns)
# cummulative distances to make angles list
cummus = np.zeros(M)
for i in range(1,M): cummus[i]=cummus[i-1]+s[i-1] 

# two lists if it starts from the middle
firstMono = M//2
sf = s[firstMono:]
sr = np.flip(s[:firstMono])
cummusf = np.zeros(len(sf)+1)
cummusr = np.zeros(len(sr)+1)
for i in range(1,len(cummusf)): cummusf[i] = cummusf[i-1] + sf[i-1]
for i in range(1,len(cummusr)): cummusr[i] = cummusr[i-1] + sr[i-1]
for i in range(M//2-1): print(sf[i],cummusf[i],sr[i],cummusr[i])
## cummulative distances to make angles list
#cummusf
#for i in range(M//2-1):
# cummusf = 

# helical coordinates lists
#angleList = cummus/rh
# starting from the middle there must be a reverse and a forward list
angleListf = cummusf/rh
angleListr = cummusr/rh
xListf = rh*np.cos(angleListf)
xListr = rh*np.cos(angleListr)
yListf = rh*np.sin(angleListf)
yListr = -rh*np.sin(angleListr)
zListf = nTurns*pitch*cummusf/sums
zListr = -nTurns*pitch*cummusr/sums

# newdmono is the list of new positions of the M monomers. The position of the first monomer must be given.
#newdmono = [np.array(firstMono)]
#for i in range(M-1):
#    newdmono.append(newdmono[i]+s[i]*r)
#newdmono = np.array(newdmono)
#print(newdmono)
# two lists if it starts from the middle
firstMonoPos = dmono[firstMono,:] # to start from the middle
newdmonof = np.zeros((len(lmonof),3))
newdmonor = np.zeros((len(lmonor),3))
helixCenter = firstMonoPos - [rh,0,0]
# puts the helix center at the origin
data = data - helixCenter
dmono = dmono - helixCenter
firstMonoPos = [rh,0,0]
newdmonof[0,:]=np.array(firstMonoPos)
newdmonor[0,:]=np.array(firstMonoPos)
for i in range(1,M//2+1):
     #print(i)
     if i<=len(sf):
         newdmonof[i,:] = [xListf[i],yListf[i],zListf[i]]
         #print(lmonof[i],newdmonof[i],angleListf[i])
     if i<=len(sr):
         newdmonor[i,:] = [xListr[i],yListr[i],zListr[i]]
         #print(lmonor[i],newdmonor[i],angleListr[i])
print(len(newdmonor),len(newdmonof))

#for i in range(M//2-1): print(lmonof[i],lmonor[i],newdmonof[i],newdmonor[i])


# slowly pulls each monomer in the list at a time.
block = 0
partialBlocks = 100
tetherparticles = []
tetherpositions = []
for part in range(M//2+3):
    print("Starting part number %d" % (part))

    a = Simulation(thermostat=0.01, temperature=300*units.kelvin)  # timestep not necessary for variableLangevin
    a.setup(platform="cuda", integrator="variableLangevin", errorTol=0.001)
    a.load(data)#, center=False)  # loads a polymer, puts a center of mass at zero
    a.saveFolder(outfolder)
 
    a.addHarmonicPolymerBonds(wiggleDist=1., bondLength=10.)
    # Bond distance will fluctuate +- 1. on average
    a.addPolynomialRepulsiveForce(trunc=5, radiusMult=10.5)
    # this accounts for topoisomeraseII allowing crossing chains occasionally (if trunc=5)
    a.addGrosbergStiffness(k=1)
  
    # fix first monomer
    if part==0: 
        tetherparticles.append(lmono[firstMono])
        tetherpositions.append(newdmonof[0,:])
    if part>0 and part<len(lmonof):
        tetherparticles.append(lmonof[part])
        tetherpositions.append(newdmonof[part,:])
    if part>0 and part<len(lmonor):
        tetherparticles.append(lmonor[part])
        tetherpositions.append(newdmonor[part,:])
    print("Fixing particles:",tetherparticles)
    print("to positions:",tetherpositions)
    a.tetherParticles(particles=tetherparticles, k=100 , positions=tetherpositions)

    if part==M/2+1: partialBlocks=50000
    if part==M/2+2:
        partialBlocks=1000        
        rcyl = 2*rh
        bottomcyl = zListr[-1]-rh
        topcyl = zListf[-1]+rh
        a.addCylindricalConfinement(r=rcyl, bottom=bottomcyl, top=topcyl, k=0.1)
        print("Cylindrical confinement has radius, top and bottom equal to %f.2, %f.2, %f.2." % (rcyl,bottomcyl,topcyl))

    # bind monomers already bound by condensins, but stop loop extrusion
    smcBondLength = 5.*a.length_scale
    wiggleDist = 2.
    k = a.kbondScalingFactor / (wiggleDist ** 2)
    for j in range(0,nLoopsII): a.forceDict["HarmonicBondForce"].addBond(int(condensinIIdata[j,0]), int(condensinIIdata[j,1]), smcBondLength, k)

    data=a.getData()
    print("data",data[0][:])
    a.step = block
    if block==0: a.save()

    # -----------Running a simulation ---------
    a.energyMinimization(stepsPerIteration=1000)
    for b in range(partialBlocks):  # Do blocks
        a.doBlock(100)  # Of 100 timesteps each
        if (a.step%100)==0: a.save()  # and save data every block

        data = a.getData()
        block = a.step

    a.save()

