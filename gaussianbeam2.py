import numpy as np                # trig functs
import scipy.constants as const   # to get physical constants
import matplotlib.pyplot as plt   # to plot
import matplotlib.patches as pat  # to draw rectangles
from matplotlib.patches import Polygon # to draw triangles

# setup parameters
numpts = 10000 # number of points to plot
displacement = np.zeros(numpts) # initialize net displacement from antenna towards plasma in mm
q = np.zeros(numpts, dtype = complex) # initialize complex beam parameter as a complex array
w = np.zeros(numpts) # initialize beam radius
freq = 94  # GHz
f   = freq*1e9
c   = const.c
Lambda  = (c/f)*1000   
x0 = 0     
xMirror = 2300   # location of flat mirror 
nhdpe = 1.52     # refractive index at 100+ GHz
w_0 = 16*0.5     # width at the antenna mouth | golsmith 169: w/a = 0.64, a = 25 | double check this number experimentally 
# TODO:
# R_h = 
# a_a = 
# w_0 = Lambda*R_h/(0.644*pi*a_a)
# RcAtMouth = R_h
RcAtMouth = 115 # 74.6125 # 115 # 47.86127
w0 = w_0/np.sqrt(1+ (((np.pi*(w_0**2))/(Lambda*RcAtMouth))**2 )) # from Goldsmitdh pg 167 beam focal radius 
print("w0: ", w0)
g_41_0  = (np.pi*(w0**2))/(Lambda)  # rayleigh length = 36.783 mm
                                    # print('Rayleigh length: ', g_41_0, 'mm')
# x_waist = RcAtMouth/(1+((Lambda*RcAtMouth)/(np.pi*(0.644*12.5)**2))**2) # print('x_waist: ', x_waist) = 31.532 mm
waveguide_length = 62 
xAn = x0 + waveguide_length # - x_waist
alfa0 = np.arctan(Lambda/(np.pi*w0)) # print('Initial divAngle: ',np.degrees(alfa0), '[deg]')
maxdisplacement =  4910 # actual is 3688, plotted longer so that you can see more of the beam path
delta_x = np.abs(maxdisplacement/numpts)

# lens parameters
# Lens 1 
thick1  = 31.2      # thickness of lens 1   
xL1     = 368 + 62  # had to add 62 to get distance to 0
RL1     = 2000      # spherical front
nL1     = nhdpe

# Lens 2 
thick2  = 60.82     # thickness of lens 2
xL2     = 1553 + 62 # 1662.5
RL2     = 1500      # radius of curvature
nL2     = nhdpe

# Lens 3 
thick3  = 51.87     # thickness of lens 3
xL3     = 1836 + 62   
RL3     = 2800
nL3     = nhdpe

# distances between elements (specifically uses the centers of lenses, not the edges)
distoL1 = np.abs(xL1-xAn) 
distoL2 = np.abs(xL2-xL1) 
distoL3 = np.abs(xL3-xL2) 

#Thick lens with first surface R1 - second surface R2, thickness (d) of material having index n2 embedded in material of index n1
f1 = RL1/(2*(nhdpe-1))
f2 = RL2/(2*(nhdpe-1))
f3 = RL3/(2*(nhdpe-1))

A1 = 1+(nhdpe-1)*thick1/(nhdpe*RL1)
A2 = 1+(nhdpe-1)*thick2/(nhdpe*RL2) 
A3 = 1+(nhdpe-1)*thick3/(nhdpe*RL3) 

B1 = thick1/nhdpe 
B2 = thick2/nhdpe 
B3 = thick3/nhdpe 

C1 = -1/f1-(nhdpe-1)**2*thick1/(nhdpe*RL1**2) 
C2 = -1/f2-(nhdpe-1)**2*thick2/(nhdpe*RL2**2) 
C3 = -1/f3-(nhdpe-1)**2*thick3/(nhdpe*RL3**2) 

D1 = 1+(1-nhdpe)*thick1/(nhdpe*RL1)
D2 = 1+(1-nhdpe)*thick2/(nhdpe*RL2)
D3 = 1+(1-nhdpe)*thick3/(nhdpe*RL3)

# start q parameter calculations
q0imag  = np.pi*(w0**2)/Lambda # q0real = R = infinity
q0      = q0imag*1j 

# p 43 Goldsmidth: q below uses thin lens approx for ABCD A=1 B=0, C= -1/f D=1
qatL1   = q0 + distoL1 

# use eqn 3.20
qaftL1  = (qatL1*A1+B1)/(qatL1*C1+D1) 
# move towards lens 2
qatL2   = qaftL1 + distoL2
qaftL2  = (qatL2*A2+B2)/(qatL2*C2+D2)
# move towards lens 3
qatL3   = qaftL2 + distoL3
qaftL3  = (qatL3*A3+B3)/(qatL3*C3+D3)

# determine values of q and w at all points during propogation through the optical system
for x in range(0,numpts-1):
    if displacement[x] <= distoL1 - thick1/2: # before lens 1
        q[x] = q0 + displacement[x]
        w[x] = np.sqrt(Lambda/(np.pi*np.imag(-1/q[x]))) # eqn 3.22b
        displacement[x+1] = displacement[x] + delta_x # increment the displacement towards the plasma
        # print('displacement[x+1] 1', displacement[x+1])
        # print('w[x] 1: ', w[x])
    elif displacement[x] > distoL1 + thick1/2 and displacement[x] <= distoL1 + distoL2 - thick2/2: # between lenses 1 and 2
        distaftL1 = displacement[x] - (distoL1 + thick1/2)
        q[x] = qaftL1 + distaftL1 # uses eqn 3.28 # fig 3
        w[x] = np.sqrt(Lambda/(np.pi*np.imag(-1/q[x]))) 
        displacement[x+1] = displacement[x] + delta_x
        # print('w[x] 1-2: ', w[x])
    
    elif displacement[x] > distoL1 + distoL2 + thick2/2 and displacement[x] <= distoL1 + distoL2 + distoL3 - thick3/2: # between lenses 2 and 3
        distaftL2 = displacement[x] - (distoL1 + distoL2 + thick2/2)
        q[x] = qaftL2 + distaftL2
        w[x] = np.sqrt(Lambda/(np.pi*np.imag(-1/q[x]))) 
        displacement[x+1] = displacement[x] + delta_x
        # print('w[x] 2-3', w[x])
    
    elif displacement[x] > distoL1 + distoL2 + distoL3 + thick3/2: # after lens 3
        distaftL3 = displacement[x] - (distoL1 + distoL2 + distoL3 + thick3/2)
        q[x] = qaftL3 + distaftL3
        w[x] = np.sqrt(Lambda/(np.pi*np.imag(-1/q[x])))
        displacement[x+1] = displacement[x] + delta_x
        # print('q[3]: ', q[x])
    
    else: # our model does not predict the q parameter and beam radius inside the lenses
        q[x] = q[x-1]
        w[x] = 0 # set equal to 0 and remove these data points from array later
        displacement[x+1] = displacement[x] + delta_x

# remove data points whose values are 0 (those inside the lenses)
displacement = displacement[w > 0]
w = w[w > 0]

# make array for position instead of displacement
position = np.zeros(displacement.size)
for x in range(0,displacement.size): # copy the old arrays over to the new arrays, shifted by position of waveguide relative to tokamak center
    position[x] = displacement[x] + xAn

# plotting parameters
z0 = 0
CL1 = [xL1+RL1-thick1/2,z0]
RL1b =  RL1
CL1b = [xL1-RL1b+thick1/2,z0]
CL2 = [xL2+RL2-thick2/2,z0]
RL2b = RL2
CL2b = [xL2-RL2b+thick2/2,z0]
CL3 = [xL3+RL3-thick3/2,z0]
RL3b = RL3
CL3b = [xL3-RL3b+thick3/2,z0]
mirrorHeight = 615
c_mirrorHeight = 203.2
ax = plt.gca()

#plot waveguides
wgLen = 58
gv = 10
wgRad = 27.0/2.0
zAn = 0 

# plot the source
whorn = 62
hhorn = 25*0.5
rec     = pat.Rectangle((x0-wgLen,z0-wgRad), wgLen, 2*wgRad, 
                        fill = False, color ='gold', zorder = 5)
ax.add_patch(rec)

pts = np.array([[0,0], [whorn, hhorn], [whorn,-hhorn]])
p = Polygon(pts, closed = False, color ='gold', zorder = 5)
ax.add_patch(p)

# plot the detector
xdet = 3816 # 3741  # 4085
print('xdet: ',xdet)
rec     = pat.Rectangle((xdet+whorn,z0-wgRad), wgLen, 2*wgRad, 
                        fill = False, color ='gold', zorder = 5)
ax.add_patch(rec)
pts = np.array([[xdet+whorn,0], [xdet, hhorn], [xdet,-hhorn]])
p = Polygon(pts, closed = False, color = 'gold', zorder = 5)
ax.add_patch(p)

# plot the lenses
#spherical surface lens 1
anglef = (np.linspace(8,-8)+180)*np.pi/180
xf = RL1*np.cos(anglef)+CL1[0]
zf = RL1*np.sin(anglef)+CL1[1]
plt.plot(xf,zf, '-', color = 'gold', zorder = 5)
#spherical back surface lens 1
anglef = (np.linspace(8,-8))*np.pi/180
xf = RL1b*np.cos(anglef)+CL1b[0]
zf = RL1b*np.sin(anglef)+CL1b[1]
plt.plot(xf,zf, '-', color = 'gold', zorder = 5)
# spherical surface lens 2
anglef = (np.linspace(12,-12)+180)*np.pi/180
xf = RL2*np.cos(anglef)+CL2[0]
zf = RL2*np.sin(anglef)+CL2[1]
plt.plot(xf,zf, '-', color = 'gold', zorder = 5)
#spherical back surface lens 2
anglef = (np.linspace(12,-12))*np.pi/180
xf = RL2b*np.cos(anglef)+CL2b[0]
zf = RL2b*np.sin(anglef)+CL2b[1]
plt.plot(xf,zf, '-', color = 'gold', zorder = 5)
#spherical surface lens 3
anglef = (np.linspace(8,-8)+180)*np.pi/180
xf = RL3*np.cos(anglef)+CL3[0]
zf = RL3*np.sin(anglef)+CL3[1]
plt.plot(xf,zf, '-', color = 'gold', zorder = 5)
#spherical back surface lens 3
anglef = (np.linspace(8,-8))*np.pi/180
xf = RL3b*np.cos(anglef)+CL3b[0]
zf = RL3b*np.sin(anglef)+CL3b[1]
plt.plot(xf,zf, '-', color = 'gold', zorder = 5)

# plot flat mirror
rec = pat.Rectangle((xMirror, z0 - (mirrorHeight/2)), 10, mirrorHeight, #40
                                fill = True, color ='hotpink')
ax.add_patch(rec) 
rec = pat.Rectangle((xMirror, z0 - (c_mirrorHeight/2)), 3.175, c_mirrorHeight, 
                                fill = True, color ='skyblue')
ax.add_patch(rec) 

# create plot
plt.plot(position, w, 'o', color='orange', markersize = 1) # plot the data
plt.plot(position, -w, 'o', color = 'orange', markersize = 1)
plt.hlines(0, whorn, xdet, color = 'silver', ls = '-.', zorder = 0)
plt.title('Beam Radius vs. Position')
plt.xlabel('Position from source [mm]') # label x-axis
plt.ylabel('1/e Electric Field Beam Radius [mm]') # label y-axis

# dataset
x_exp = [380, 380, 380, 1570, 1570, 1570, 1860, 1860, 1860]
y_exp = [ 38,  47,  33,  186,  186,  195,  174,  174,  182]
x_avg = [380, 1570, 1860, 3741, 3741]
y_avg = [ 39,  189,  177,   28,   45]
plt.plot(x_exp, y_exp,'o', color = 'skyblue', markersize = 3, zorder = 5)
plt.plot(x_avg, y_avg,'o', color = 'tomato', markersize = 3, zorder = 7)
plt.errorbar(x_exp, y_exp, xerr = 5, yerr = 5, ls = 'none', ecolor = 'paleturquoise', capsize = 2, elinewidth = 0.5, barsabove = False)
plt.errorbar(x_avg, y_avg, xerr = 5, yerr = 5, ls = 'none', ecolor = 'firebrick',     capsize = 2, elinewidth = 0.5, barsabove = False, zorder = 6)

plt.show()

# plotting what occurs for the 8" mirror situation
# same lens 1, 2 equations
# just draw a line connecting back of lens 2 to len 3?
# after lens 3, i need an equation, have data points to fit it
# so the beam comes from lens 2, and gets cut off by 48.26 mm on the radius
# equivalent to beam being blocked by blocking material with an 8" diameter hole cut out, yes?
# so that means that the gaussian beam get truncated - see on the graph
# to do: try to find equation for beam getting truncated - or try to create own