import numpy as np                # trig functs
import scipy.constants as const   # to get physical constants
import matplotlib.pyplot as plt   # to plot
import matplotlib.patches as pat  # to draw rectangles
from matplotlib.patches import Polygon # to draw triangles

'''
TODO:
- 
'''
lines = []
neg_lines = []
freq_list = []
w_master_list = []
position_master_list = []
color_master_list = []
x0 = 0    #set 0 as plasma center
a = 570
numpts = 100000 # number of points to plot
c   = const.c
xWindow = 2390      # location of window 
nmetal = 2.75681
w_0 = 8   # width at the plasma
distance = 27788 #27823 # 35427 # 35455 #22027
xpos = []
point_pos = 27764 # distance
#   fs = RLs/(2*(nmetal-1))
thicknesses = np.array([5, 5, 5, 5,     5,     5,     5,     5])
xLs = np.array([1501, 3900, 6150, 11027, 20430, 26858, 27155, 27655])
#xLs = np.array([1501, 3900, 6150, 7927, 11027, 20430, 26858, 27155, 27655])

'''
disbt0 = xLs[1] - xLs[0] # 3900-1501
findRL0 = (disbt0-(2500/(2*(nmetal-1))))*(2*(nmetal-1))
print(f'RL2 = {findRL0}')
dis0 = 2500/(2*(nmetal-1)) + findRL0/(2*(nmetal-1))
print(f'dis0: {disbt0, dis0}')
'''

disbt1 = xLs[2] - xLs[1] 
findRL1 = (disbt1-(2500/(2*(nmetal-1))))*(2*(nmetal-1))
print(f'RL1 = {findRL1}')    # 5406
dis1 = 2500/(2*(nmetal-1)) + findRL1/(2*(nmetal-1))
print(f'dis1: {disbt1, dis1}')

disbt2 = xLs[3] - xLs[2] 
findRL2 = (disbt2-(findRL1/(2*(nmetal-1))))*(2*(nmetal-1))
print(f'RL2 = {findRL2}')    # 838
dis2 = findRL1/(2*(nmetal-1)) + findRL2/(2*(nmetal-1))
print(f'dis2: {disbt2, dis2}')

disbt3 = xLs[4] - xLs[3] 
findRL3 = (disbt3-(findRL2/(2*(nmetal-1))))*(2*(nmetal-1))
print(f'RL3 = {findRL3}')    # 10054
dis3 = findRL3/(2*(nmetal-1)) + findRL2/(2*(nmetal-1))
print(f'dis3: {disbt3, dis3}')

disbt4 = xLs[5] - xLs[4] 
findRL4 = (disbt4-(findRL3/(2*(nmetal-1))))*(2*(nmetal-1))
print(f'RL4 = {findRL4}')    # 
dis4 = findRL4/(2*(nmetal-1)) + findRL3/(2*(nmetal-1))
print(f'dis4: {disbt4, dis4}')

disbt5 = xLs[6] - xLs[5] 
findRL5 = (disbt5-(findRL4/(2*(nmetal-1))))*(2*(nmetal-1))
print(f'RL5 = {findRL5}')    # 
dis5 = findRL5/(2*(nmetal-1)) + findRL4/(2*(nmetal-1))
print(f'dis5: {disbt5, dis5}')

disbt6 = xLs[7] - xLs[6] 
findRL6 = (disbt6-(findRL5/(2*(nmetal-1))))*(2*(nmetal-1))
print(f'RL6 = {findRL6}')    # 
dis6 = findRL6/(2*(nmetal-1)) + findRL5/(2*(nmetal-1))
print(f'dis6: {disbt6, dis6}')

'''
disbt7 = xLs[8] - xLs[7] 
findRL7 = (disbt7-(findRL6/(2*(nmetal-1))))*(2*(nmetal-1))
print(f'RL7 = {findRL7}')    # 10054
dis7 = findRL7/(2*(nmetal-1)) + findRL6/(2*(nmetal-1))
print(f'dis7: {disbt7, dis7}')
'''

disbt7 = point_pos - xLs[7] 
findRL7 = (disbt7-(findRL6/(2*(nmetal-1))))*(2*(nmetal-1))
print(f'RL7 = {findRL7}')    # 10054
dis7 = findRL7/(2*(nmetal-1)) + findRL6/(2*(nmetal-1))
print(f'dis7: {disbt7, dis7}')

'''
disbt8 = point_pos - xLs[8] 
findRL8 = (disbt8-(findRL7/(2*(nmetal-1))))*(2*(nmetal-1))
print(f'RL8 = {findRL8}')    # 10054
dis8 = findRL8/(2*(nmetal-1)) + findRL7/(2*(nmetal-1))
print(f'dis8: {disbt8, dis8}')
'''

#RLs = np.array([2500, 3000, 5000, 7000, 7000, 8000, 7000, 7461, 450])
RLs = np.array([2500, findRL1, 15000, findRL3, 20000, 10000, 12000, 400])
#                             10000   21308         
# basically same waist at detector, but feasable spot to put the detector 
# is within 3mm; all slope is steep, so v hard not to clip

f0 = RLs[0]/(2*(nmetal-1))
f1 = RLs[1]/(2*(nmetal-1))
f2 = RLs[2]/(2*(nmetal-1))
f3 = RLs[3]/(2*(nmetal-1))
f4 = RLs[4]/(2*(nmetal-1))
f5 = RLs[5]/(2*(nmetal-1))
f6 = RLs[6]/(2*(nmetal-1))
f7 = RLs[7]/(2*(nmetal-1))
# f8 = RLs[8]/(2*(nmetal-1))

findis1 = f0 + f1 + 1501    # 3859.8208172767686 for rl2 = 5788
findis2 = f1 + f2 + findis1
findis3 = f2 + f3 + findis2 
findis4 = f3 + f4 + findis3 
findis5 = f4 + f5 + findis4 
findis6 = f5 + f6 + findis5 
findis7 = f6 + f7 + findis6
#findis8 = f7 + f8 + findis7 

# 26838.36; 26689.7
number_lenses = len(xLs)


def beam_prop(freq, color_line, xAn):
    q = np.zeros(numpts, dtype = complex) # initialize complex beam parameter as a complex array
    w = np.zeros(numpts)
    displacement = np.zeros(numpts) # initialize net displacement from antenna towards plasma in mm
    maxdisplacement = distance + a
    delta_x = np.abs(maxdisplacement/numpts)
    f   = freq*1e9
    Lambda  = (c/f)*1000   # units in mm 
    xAn = xAn
    print('freq:', freq)

    # sanity check
    if len(RLs) == len(xLs) == len(thicknesses):
        pass
    else:
        raise Exception(f'you done goofed up kiddo: len(RLs) != len(xLs) != len(thicknesses)')
    
    
    pre_xLRH = np.array([xAn])
    xLRH = np.append(pre_xLRH, xLs[:number_lenses-1])
    B = np.ones(number_lenses)
    before_distances = []
    after_distances = []

    q0imag  = np.pi*(w_0**2)/Lambda # q0real = R = infinity
    q0      = q0imag*1j
    q_at_values = []
    qaft = [q0]
    qs_at = {}
    qs_aft = {}
    xpos_names = [xAn]

    # calculations
    disto = np.abs(xLs - xLRH)
    fs = RLs/(2*(nmetal-1))
    A = 1+(nmetal-1)*thicknesses/(nmetal*RLs)
    B = B*thicknesses/nmetal 
    C = -1/fs -(nmetal-1)**2*thicknesses/(nmetal*RLs**2)
    D = 1+(1-nmetal)*thicknesses/(nmetal*RLs)

    for i in range(number_lenses): 
        # q parameter calculations
        q_at_name = f'qatL{str(i)}'
        q_aft_name = f'qaftL{str(i)}'
        last_q_value = qaft[-1]  
        last_dist_value = disto[i]    # get latest disto value
        qs_at[q_at_name] = last_q_value + last_dist_value
        q_at_values.append(last_q_value + last_dist_value)
        q_after = (q_at_values[i]*A[i]+B[i])/(q_at_values[i]*C[i]+D[i])
        qs_aft[q_aft_name] = q_after
        qaft.append(q_after)    

        # distances for the next part
        up_to_sum = np.sum(disto[:i + 1])
        before_distances.append(up_to_sum - thicknesses[i]/2)
        after_distances.append(up_to_sum + thicknesses[i]/2)

    # determine values of q and w at all points during propogation through the optical system
    for x in range(0,numpts-1):
        if displacement[x] <= before_distances[0]: 
            q[x] = q0 + displacement[x]
            w[x] = np.sqrt(Lambda/(np.pi*np.imag(-1/q[x]))) # eqn 3.22b
            displacement[x+1] = displacement[x] + delta_x # increment the displacement towards the plasma
       
        elif displacement[x] > after_distances[0] and displacement[x] <= before_distances[1]: # between lenses 1 and 2 
            distaftL0 = displacement[x] - (after_distances[0])
            q[x] = qaft[1] + distaftL0 # uses eqn 3.28 # fig 3
            w[x] = np.sqrt(Lambda/(np.pi*np.imag(-1/q[x]))) 
            displacement[x+1] = displacement[x] + delta_x
       
        elif displacement[x] > after_distances[1] and displacement[x] <= before_distances[2]: # between lenses 2 and 3
            distaftL1 = displacement[x] - (after_distances[1])
            q[x] = qaft[2] + distaftL1
            w[x] = np.sqrt(Lambda/(np.pi*np.imag(-1/q[x]))) 
            displacement[x+1] = displacement[x] + delta_x
       
        elif displacement[x] > after_distances[2] and displacement [x] <= before_distances[3]: # between lenses 3 & 4
            distaftL2 = displacement[x] - (after_distances[2])
            q[x] = qaft[3] + distaftL2
            w[x] = np.sqrt(Lambda/(np.pi*np.imag(-1/q[x])))
            displacement[x+1] = displacement[x] + delta_x
       
        elif displacement[x] > after_distances[3] and displacement[x] <= before_distances[4]: # between lens 4 & 5
            distaftL3 = displacement[x] - (after_distances[3])
            q[x] = qaft[4] + distaftL3
            w[x] = np.sqrt(Lambda/(np.pi*np.imag(-1/q[x])))
            displacement[x+1] = displacement[x] + delta_x

        elif displacement[x] > after_distances[4] and displacement[x] <= before_distances[5]: # after lens 5
            distaftL4 = displacement[x] - (after_distances[4])
            q[x] = qaft[5] + distaftL4
            w[x] = np.sqrt(Lambda/(np.pi*np.imag(-1/q[x])))
            displacement[x+1] = displacement[x] + delta_x

        elif displacement[x] > after_distances[5] and displacement[x] <= before_distances[6]: # after lens 5
            distaftL5 = displacement[x] - (after_distances[5])
            q[x] = qaft[6] + distaftL5
            w[x] = np.sqrt(Lambda/(np.pi*np.imag(-1/q[x])))
            displacement[x+1] = displacement[x] + delta_x
        
        elif displacement[x] > after_distances[6] and displacement[x] <= before_distances[7]: # after lens 5
            distaftL6 = displacement[x] - (after_distances[6])
            q[x] = qaft[7] + distaftL6
            w[x] = np.sqrt(Lambda/(np.pi*np.imag(-1/q[x])))
            displacement[x+1] = displacement[x] + delta_x
        
        elif displacement[x] > after_distances[7]: #and displacement[x] <= before_distances[8]: # after lens 5
            distaftL7 = displacement[x] - (after_distances[7])
            q[x] = qaft[8] + distaftL7
            w[x] = np.sqrt(Lambda/(np.pi*np.imag(-1/q[x])))
            displacement[x+1] = displacement[x] + delta_x
            '''
        elif displacement[x] > after_distances[8]: # after lens 5
            distaftL8 = displacement[x] - (after_distances[8])
            q[x] = qaft[9] + distaftL8
            w[x] = np.sqrt(Lambda/(np.pi*np.imag(-1/q[x])))
            displacement[x+1] = displacement[x] + delta_x
            '''
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
    CLs = [xLs+RLs-thicknesses/2]  # would be [xLs+RLs-thicknesses/2, z0] but since z = 0, omitted
    RLbs =  RLs
    CLbs = [xLs-RLbs+thicknesses/2]


    windowHeight = 254
    ax = plt.gca()

    #plot waveguides
    wgLen = 46
    gv = 10
    wgRad = 27.0/2.0
    zAn = 0 

    # plot the detector
    global whorn
    ### EDIT HERE ###
    whorn = 46
    hhorn = 4.6*0.5
    global xdet
    xdet = point_pos  
    rec = pat.Rectangle((xdet+whorn,z0-wgRad), wgLen, 2*wgRad, 
                    fill = False, color ='gold')   
    ax.add_patch(rec)
    pts = np.array([[xdet+whorn,0], [xdet, hhorn], [xdet,-hhorn]])
    p = Polygon(pts, closed = False, color = 'gold') 
    ax.add_patch(p)


    # plot the lenses
    #spherical surfaces
    for i in range(number_lenses):
        anglef = (np.linspace(8,-8)+180)*np.pi/180
        xf = RLs[i]*np.cos(anglef) + CLs[0][i] # RLs[i]
        zf = RLs[i]*np.sin(anglef)    # would be + CL[1]
        plt.plot(xf, zf, '-', color = 'gold', zorder = 5)

    #spherical back surfaces
    for i in range(number_lenses):
        anglef = (np.linspace(8,-8))*np.pi/180
        xf = RLbs[i]*np.cos(anglef)+CLbs[0][i]
        zf = RLbs[i]*np.sin(anglef)    # would be + CL[1]
        plt.plot(xf, zf, '-', color = 'gold', zorder = 5)

    # plot wplasma
    circle = plt.Circle((0, 0), a, color = 'lightcoral')    # plot plasma
    ax.add_patch(circle)

    # plot window
    rec = pat.Rectangle((xWindow, z0 - (windowHeight/2)), 40, windowHeight,
                                fill = True, color ='mediumpurple')
    ax.add_patch(rec) 

    #plot ICRF antenna
    xicrf = 4705
    rec = pat.Rectangle((xicrf, z0 - (windowHeight/2)), 40, windowHeight,
                                fill = True, color ='mediumpurple')
    ax.add_patch(rec) 

    # plot width limit
    plt.axhline(y = 50.8, color = 'mediumpurple', zorder = 0) # 2" limit
    plt.axhline(y = -50.8, color = 'mediumpurple', zorder = 0)
    plt.axhline(y = 44.45, color = 'mediumpurple', zorder = 0)  # 1.75" limit for <7" beam bt L1 and window
    plt.axhline(y = -44.45, color = 'mediumpurple', zorder = 0)

    #plt.ylim(-200,200)
    freq_list.append(freq)
    # print('freq_list', freq_list)

    w_list = []
    position_list = []
    for i in range(len(w)):
        w_list.append(w[i])
    w_master_list.append(w_list)
    for i in range(len(position)):
        position_list.append(position[i])
    position_master_list.append(position_list)
    color_master_list.append

    for i in range(len(freq_list)):
        plt.plot(position, w, 'o', color = color_line, markersize = 1) #, label = freq_list[-1]) # plot the data
        plt.plot(position, -w, 'o', color = color_line, markersize = 1)
        plt.hlines(0, whorn, xdet, color = 'silver', ls = '-.', zorder = 0)
        plt.title('Beam Radius vs. Position')
        plt.xlabel('Position from source [mm]') # label x-axis
        plt.ylabel('1/e Electric Field Beam Radius [mm]') # label y-axis

    # dataset
    x_exp = [point_pos]
    y_exp = [1.25] 
    plt.plot(x_exp, y_exp,'o', color = 'skyblue', markersize = 3, zorder = 5)

    plt.ylim(-100,100)
    plt.xlim(-a, distance + a)

# 360-470; 245-345
beam_prop(245, 'red', x0 + a)    # 240
beam_prop(260, 'orange', x0 + a)    # 240
beam_prop(340, 'mediumblue', x0)    #360
beam_prop(470, 'purple', x0)    # 450

plt.show()