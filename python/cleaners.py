import numpy as np  
import math
    
def noisered(m0,m1,m,l,s):
    #m0 = moment 0 map, m1 = moment 1 map, m = max x pixels, l = max y pixels, s = strenght of cleaning algorithm
    for n in range(100):
        for i in range(0,m):
            for j in range(0,l):
                flag = False
                if m0[i,j]!= 0.0:
                    for k in range(1,s): #finds noise showiging up on mom1 map, and flags it
                        if i < (m-k) and j < (l-k) and i > k and j > k:
                            if m0[i+k,j+k] == 0.0 and m0[i-k,j+k] == 0.0 and m0[i+k,j-k] == 0.0 and m0[i-k,j-k] == 0.0:
                                flag = True
                            if m0[i+k,j] == 0.0 and m0[i-k,j] == 0.0 and m0[i,j-k] == 0.0 and m0[i,j+k] == 0.0:
                                flag = True
                if flag:
                    m0[i,j] = 0.0
                    m1[i,j] = math.nan
                    
                    
def cleancube(m0,c,m,l,v,no0=False):
    #m0 = moment 0 map, c = cube, m = max x pixels, l = max y pixels, v = nr of channels
    for i in range(0,m):
        for j in range(0,l):
            for o in range(0,v+1):
                if m0[i,j] == 0: #zeroez pixels that are zero in mom0 on cube
                    c[o,i,j] = 0.0
                if no0:
                    if c[o,i,j] < 0:
                        c[o,i,j] = 0.0
                    
                    
def cleanloop(m0,m1,c,m,l,v,vmin,vmax,s, nan = False, cube=False): 
    #m0 = moment 0 map, m1 = moment 1 map, c = cube, m = max x pixels, l = max y pixels
    #vmin = minimum velocity, vmax = maximum velocity, s = strenght of cleaning algorithm
    print('starting cleanup')
    for i in range(0, m):
        for j in range(0, l):
            
            if m0[i,j] <= 0.0: # sets bounds of galaxy from mom 0 map, anything not showing up on mom 0 will be flagged
                flag = True
            else:
                flag = False
                            
                if m1[i,j] > vmax or m1[i,j] < vmin:
                    flag = True
                    
            if flag: #removes all flagged pixels from mom0 and mom1
                m0[i,j] = 0.0
                m1[i,j] = math.nan
    
    print('cleanup finished, starting noise reduction')
    noisered(m0,m1,m,l,s) #noise reduction algorithm           
    print('noise reduction finished')
    
    if cube:
        print('cleaning cube')           
        cleancube(m0,c,m,l,v)#cube cleaning/masking algorithm
        print('cube cleaned')
    
    if nan:
        for i in range(0, m):
            for j in range(0, l):
                if m0[i,j] == 0.0:
                    m0[i,j] = math.nan
                    
def createmask(inm, m0 ,m,l,v):
    print('creating mask')
    outm = np.zeros([v+1,m,l])
    for i in range(0,m):
        for j in range(0,l):
            for k in range(0,v+1):
                if inm[k,i,j] > 0:
                    outm[k,i,j] = 1
                if m0[i,j] == 0:
                    outm[k,i,j] = 0
    print('mask made, cleaning mask')
    cleanmask(outm,m,l,v)
    print('mask cleaned')
    return outm

def cleanmask(inm,m,l,v):
     for n in range(100):
        for o in range(0,v+1):
            for i in range(0,m):
                for j in range(0,l):
                    flag = False
                    if inm[o,i,j]!= 0.0:
                        for k in range(1,4): #finds noise showiging up on mom1 map, and flags it
                            if i < (m-k) and j < (l-k) and i > k and j > k:
                                if inm[o,i+k,j+k] == 0.0 and inm[o,i-k,j+k] == 0.0 and inm[o,i+k,j-k] == 0.0 and inm[o,i-k,j-k] == 0.0:
                                    flag = True
                                if inm[o,i+k,j] == 0.0 and inm[o,i-k,j] == 0.0 and inm[o,i,j-k] == 0.0 and inm[o,i,j+k] == 0.0:
                                    flag = True
                if flag:
                    inm[o,i,j] = 0.0
    