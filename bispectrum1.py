from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.stats import norm
import matplotlib
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.mlab as mlab
from astropy.io import fits
from math import pi, sin, cos, sqrt, log, floor
from sympy.physics.wigner import gaunt
import sys
from utils import to_str
from main import gaunt_gen
from pipeline import gaunt_oj
import pipeline

l_dic ={}
gaunts=0

def gauntB(l1, l2, l3, m1, m2, m3, p1, p2, p3, alm1, alm2, alm3):
    #function to change
    b_ell = gaunt(l1,l2,l3,m1,m2,m3)*np.real(alm1[p1]*alm2[p2]*alm3[p3])
    return b_ell

# def gauntB_oj(l1, l2, l3, m1, m2, m3, p1, p2, p3, alm1, alm2, alm3):
#     l_arr=np.array([l1, l2, l3])
#     m_arr=np.array([m1, m2, m3])
#     l_idx=to_str(l_arr)
#     if l_idx not in l_dic.keys():
#         gaunt_gen(l_arr)
#         l_dic[l_idx]=1
#     b_ell=gaunt_oj(l_arr,m_arr)*np.real(alm1[p1]*alm2[p2]*alm3[p3])    
#     return b_ell

def gauntB_oj(l1, l2, l3, m1, m2, m3):
    l_arr=np.array([l1, l2, l3])
    m_arr=np.array([m1, m2, m3])
    l_idx=to_str(l_arr)
    if l_idx not in l_dic.keys():
        gaunt_gen(l_arr)
        l_dic[l_idx]=1
    b_ell=gaunt_oj(l_arr,m_arr)   
    return b_ell

def alm_position(lmax, la, ma):
    alm_position = hp.Alm.getidx(lmax, la, ma)
    return alm_position

def equilatero(i):
    global gaunts
    #### Select your .fits files (map file) #############
    mapa = hp.read_map()
    #mask = hp.read_map()
    #mapa = mapa*mask
    #######################################
     
    lmax = 300   ## l_max you need to calculate
    lmaxm = 767  ## l_max of your map

    B_ell = []
    l1_ell = []
    l2_ell = []
    l3_ell = []
          
    alm1 = hp.map2alm(mapa)
    alm2 = hp.map2alm(mapa)
    alm3 = hp.map2alm(mapa)

    l1 = 0
    l2 = 0
    l3 = 0

    for l1 in range(0, lmax+1):                   
        sum=0
        if ((l1+l2+l3)%2 == 0):
            l2=l1
            l3=l1
            l1_ell.append(l1)
            l2_ell.append(l2)
            l3_ell.append(l3)
            for m1 in range(-l1, l1+1):
                p1 = alm_position(lmaxm, l1, abs(m1))
                x=alm1[p1]
                for m2 in range(-l2, l2+1):                    
                    p2 = alm_position(lmaxm, l2, abs(m2))
                    y=alm2[p2]
                    for m3 in range(-l3, l3+1):                        
                        p3 = alm_position(lmaxm, l3, abs(m3))
                        z=alm3[p3]
                        # sum =float(sum)+ gauntB(l1, l2, l3, m1, m2, m3, p1, p2, p3, alm1, alm2, alm3)
                        if m1+m2+m3==0:
                            gaunts+=1
                            sum = float(sum)+ gauntB_oj(l1, l2, l3, m1, m2, m3)*np.real(x*y*z)
            B_ell.append(sum)
            pipeline.gaunt_dic.clear()
            l_dic.clear()
    bl = np.array(B_ell)
    l1 = np.array(l1_ell)
    l2 = np.array(l2_ell)
    l3 = np.array(l3_ell)
    bl = bl.astype(float)
    l1 = l1.astype(float)
    l2 = l2.astype(float)
    l3 = l3.astype(float)
    result = np.vstack([bl, l1, l2, l3])

    np.savetxt('name_of_your_file_equilateral.csv', result.T, delimiter = ',', header = 'bl, l1, l2, l3')

def equisize(i):
    global gaunts
    
    #### Select your .fits files (map file) #############
    mapa = hp.read_map()
    #mask = hp.read_map()
    #mapa = mapa*mask
    #######################################

    lmax = 300     ## l_max you need to calculate
    lmaxm = 767    ## l_max of your map

    B_ell = []
    l1_ell = []
    l2_ell = []
    l3_ell = []
          
    alm1 = hp.map2alm(mapa)
    alm2 = hp.map2alm(mapa)
    alm3 = hp.map2alm(mapa)
    
    l0 = lmax
    l1 = 0
    l2 = 0
    l3 = 0
        

        # print("Entrei nos fors")
    for l1 in range(0, lmax+1):
        for l2 in range(0, lmax+1):
            for l3 in range(0, lmax+1):
                if((l1<=l2<=l3) and (abs(l1-l2)<=l3<=abs(l1+l2)) and ((l1+l2+l3)==l0)):
                    l1_ell.append(l1)
                    l2_ell.append(l2)
                    l3_ell.append(l3)

                        #call here                   
                    sum=0
                    for m1 in range(-l1, l1+1):
                        p1 = alm_position(lmaxm, l1, abs(m1))
                        x=alm1[p1]
                        for m2 in range(-l2, l2+1):                    
                            p2 = alm_position(lmaxm, l2, abs(m2))
                            y=alm2[p2]
                            for m3 in range(-l3, l3+1):                        
                                p3 = alm_position(lmaxm, l3, abs(m3))
                                z=alm3[p3]
                                # sum =float(sum)+ gauntB(l1, l2, l3, m1, m2, m3, p1, p2, p3, alm1, alm2, alm3)
                                if m1+m2+m3==0:
                                    gaunts+=1
                                    sum = float(sum)+ gauntB_oj(l1, l2, l3, m1, m2, m3)*np.real(x*y*z)
                                    

                        
                    B_ell.append(sum)
                    pipeline.gaunt_dic.clear()
                    l_dic.clear()
    bl = np.array(B_ell)
    l1 = np.array(l1_ell)
    l2 = np.array(l2_ell)
    l3 = np.array(l3_ell)
    bl = bl.astype(float)
    l1 = l1.astype(float)
    l2 = l2.astype(float)
    l3 = l3.astype(float)
    result = np.vstack([bl, l1, l2, l3])

    np.savetxt('name_of_your_file_equisize.csv', result.T, delimiter = ',', header = 'bl, l1, l2, l3')

def isosceles(i):
    global gaunts
    
    #### Select your .fits files (map file) #############
    mapa = hp.read_map()
    #mask = hp.read_map()
    #mapa = mapa*mask
    #######################################
    
    
    lmax = 150     ## l_max you need to calculate
    lmaxm = 767    ## l_max of your map

    B_ell = []
    l1_ell = []
    l2_ell = []
    l3_ell = []

          
    alm1 = hp.map2alm(mapa)
    alm2 = hp.map2alm(mapa)
    alm3 = hp.map2alm(mapa)
    
    l0 = lmax
    l1 = 0
    l2 = 0
    l3 = 0
        

        # print("Entrei nos fors")
    for l1 in range(0, lmax+1):
        for l2 in range(0, lmax+1):
            for l3 in range(0, lmax+1):
                if((l1<=l2<=l3) and (abs(l1-l2)<=l3<=abs(l1+l2)) and (l1==l2) and ((l1+l2+l3)%2 == 0)):
                    l1_ell.append(l1)
                    l2_ell.append(l2)
                    l3_ell.append(l3)

                        #call here                   
                    sum=0
                    for m1 in range(-l1, l1+1):
                        p1 = alm_position(lmaxm, l1, abs(m1))
                        x=alm1[p1]
                        for m2 in range(-l2, l2+1):                    
                            p2 = alm_position(lmaxm, l2, abs(m2))
                            y=alm2[p2]
                            for m3 in range(-l3, l3+1):                        
                                p3 = alm_position(lmaxm, l3, abs(m3))
                                z=alm3[p3]
                                # sum =float(sum)+ gauntB(l1, l2, l3, m1, m2, m3, p1, p2, p3, alm1, alm2, alm3)
                                if m1+m2+m3==0:
                                    gaunts+=1
                                    sum = float(sum)+ gauntB_oj(l1, l2, l3, m1, m2, m3)*np.real(x*y*z)
                                    

                        
                    B_ell.append(sum)
                    pipeline.gaunt_dic.clear()
                    l_dic.clear()
    bl = np.array(B_ell)
    l1 = np.array(l1_ell)
    l2 = np.array(l2_ell)
    l3 = np.array(l3_ell)
    bl = bl.astype(float)
    l1 = l1.astype(float)
    l2 = l2.astype(float)
    l3 = l3.astype(float)
    result = np.vstack([bl, l1, l2, l3])

    np.savetxt('name_of_your_file_isosceles.csv', result.T, delimiter = ',', header = 'bl, l1, l2, l3')
