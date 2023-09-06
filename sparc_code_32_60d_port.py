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
x0 = 0    # set 0 as edge of horn
a = 570
numpts = 100000 # number of points to plot
c   = const.c
xWindow = 2390      # location of window 
nmetal = 2.75681
distance = 22686 #21756 #23465
xpos = []
point_pos = distance
#   fs = RLs/(2*(nmetal-1))

thicknesses = np.array([5, 5, 5, 5])
sec_m = 18780
xLs = np.array([distance-(1509+295+6372+10630+2740+1000), distance-(1509+295+6372+10630), distance - (1509+295+6372), distance - 1509])
RLs = np.array([500, 10000, 15000, 3900])

print(RLs/(2*(nmetal-1)))
print(500/(2*(nmetal-1)))

'''new distances'''
# RLs = np.array([500, 10000, 15000, 3900]) # 10.5, 15.5, sent to nathan
# RLs = np.array([500, 10000, 15000, 4100]) #11.4, 14.2


'''old distances'''
# RLs = np.array([500, 10000, 15000, 3900]) #9.6, 17
# RLs = np.array([500, 10000, 13000, 4000]) #12, 15
# RLs = np.array([450, 10000, 13000, 2800]) #13, 16 | what i sent to nathan
# RLs = np.array([450, 10000, 12000, 2800]) #13, 16
# RLs = np.array([450, 13000, 13000, 2200]) #11, 18
# RLs = np.array([450, 12000, 13000, 2400]) 



number_lenses = len(xLs)


def beam_prop(freq, color_line, xAn, w_0):
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
       
        elif displacement[x] > after_distances[3]: #and displacement[x] <= before_distances[4]: # between lens 4 & 5
            distaftL3 = displacement[x] - (after_distances[3])
            q[x] = qaft[4] + distaftL3
            w[x] = np.sqrt(Lambda/(np.pi*np.imag(-1/q[x])))
            displacement[x+1] = displacement[x] + delta_x
            '''
        elif displacement[x] > after_distances[4] # and displacement[x] <= before_distances[5]: # between lens 5 & 6
            distaftL4 = displacement[x] - (after_distances[4])
            q[x] = qaft[5] + distaftL4
            w[x] = np.sqrt(Lambda/(np.pi*np.imag(-1/q[x])))
            displacement[x+1] = displacement[x] + delta_x
            
        elif displacement[x] > after_distances[5]: and displacement[x] <= before_distances[6]: # between lens 6 & 7
            distaftL5 = displacement[x] - (after_distances[5])
            q[x] = qaft[6] + distaftL5
            w[x] = np.sqrt(Lambda/(np.pi*np.imag(-1/q[x])))
            displacement[x+1] = displacement[x] + delta_x
            
        elif displacement[x] > after_distances[6] and displacement[x] <= before_distances[7]: # between lens 7 & 8
            distaftL6 = displacement[x] - (after_distances[6])
            q[x] = qaft[7] + distaftL6
            w[x] = np.sqrt(Lambda/(np.pi*np.imag(-1/q[x])))
            displacement[x+1] = displacement[x] + delta_x
            
        elif displacement[x] > after_distances[7]: # and displacement[x] <= before_distances[8]: # after lens 8
            distaftL7 = displacement[x] - (after_distances[7])
            q[x] = qaft[8] + distaftL7
            w[x] = np.sqrt(Lambda/(np.pi*np.imag(-1/q[x])))
            displacement[x+1] = displacement[x] + delta_x
            
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
    xdet = -46 # x0-wgLen
    rec = pat.Rectangle((2*xdet,z0-wgRad), wgLen, 2*wgRad, 
                            fill = False, color ='gold')
    ax.add_patch(rec)
    pts = np.array([[xdet,0], [0, 4.6*0.5], [0,-4.6*0.5]])
    p = Polygon(pts, closed = False, color ='gold')
    ax.add_patch(p)  
    circle = plt.Circle((distance - 0, 0), a, color='lightcoral') # plot plasma
    ax.add_patch(circle)


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

    # plot window
    rec = pat.Rectangle((distance - xWindow, z0 - (windowHeight/2)), 40, windowHeight,
                                fill = True, color ='mediumpurple')
    ax.add_patch(rec) 

    #plot ICRF antenna
    xicrf = 4705
    rec = pat.Rectangle((distance - xicrf, z0 - (windowHeight/2)), 40, windowHeight,
                                fill = True, color ='mediumpurple')
    ax.add_patch(rec) 

    # plot lines in plasma
    rec = pat.Rectangle((distance, z0 - (windowHeight/2)), 4, windowHeight,
                                fill = True, color ='crimson')
    ax.add_patch(rec) 

    rec = pat.Rectangle((distance - a, z0 - (windowHeight/2)), 4, windowHeight,
                                fill = True, color ='crimson')
    ax.add_patch(rec) 

    # plot width limit
    plt.axhline(y = 50.8, color = 'mediumpurple', zorder = 0) # 2" limit
    plt.axhline(y = -50.8, color = 'mediumpurple', zorder = 0)
    plt.axhline(y = 44.45, color = 'mediumpurple', zorder = 0)  # 1.75" limit for <7" beam bt L1 and window
    plt.axhline(y = -44.45, color = 'mediumpurple', zorder = 0)

    plt.ylim(-100,100)
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
    x_exp = [point_pos, point_pos - a]
    y_exp = [15,         8] 
    plt.plot(x_exp, y_exp,'o', color = 'skyblue', markersize = 3, zorder = 5)

    #plt.ylim(-100,100)
    plt.xlim(-a, distance + a)

# 360-470; 245-345
beam_prop(260, 'red', x0, w_0 = 1.9)    # 240
beam_prop(340, 'orange', x0, w_0 = 1.5)    # 240
beam_prop(340, 'mediumblue', x0, w_0 = 1.9)    #360
beam_prop(450, 'purple', x0, w_0 = 1.5)    # 450
beam_prop(245, 'hotpink', x0, w_0 = 1.9)

plt.show()