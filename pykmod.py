#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 20:54:57 2020

@author: jared
"""

import numpy as np
from scipy.integrate import simpson
from scipy.interpolate import RegularGridInterpolator
from scipy.interpolate import interp1d
import sys
import os.path

#kpath='/home/jared/fl/fl_0.79.47/lib/'

availteff=np.append(np.append(np.asarray([3500+(250*i) for i in range(35)]),np.asarray([12500+(500*i) for i in range(16)])),np.asarray([21000+(1000*i) for i in range(10)]))
availlogg=np.asarray([0.5*i for i in range(11)])
availmetal=np.asarray([-5,-4.5,-4,-3.5,-3,-2.8,-2.5,-2.3,-2.0,-1.8,-1.5,-1.3,-1.0,-0.8,-0.5,-0.3,0,0.3,0.5,0.8,1,1.5])
#availmetal=np.array([-2.5+0.5*i for i in range(10)])

def rd_kmod(teff,logg,metal):                                                                 
    if metal >= 0:                                               
        s2='p'                                                  
    else:                                                        
        s2='m'                                                  
                                                                 
    s3 = str(metal).replace('.','').replace('-','')     
         
    path = os.path.abspath(os.path.dirname(__file__))
    kpath= path + '/modelatmospheres/'                                                         
    #filename = kpath+'a'+s2+s3+'k2odfnew.dat'                    
    filename = kpath+'m'+s2+s3+'cp00op00.mod'  
                                                           
    text = 'texto'                                               
    teffstring="{:7.0f}".format(teff)                            
    loggstring="{:8.5f}".format(logg)                            
                                                                 
    file = open(filename, 'r')                                   
    text = file.readlines()                                      
                                                                 
    for i in range(len(text)):                                   
        if teffstring in text[i] and loggstring in text[i]:      
            header = text[i:i+23]                                
            model = np.zeros((72,10))                            
            for j in range(72):                                  
                model[j] = text[i+23+j].split()                  
            footer = text[i+95:i+97]                             
            break                                                
    file.close()                                                 
    return (model,header,footer)

def pykmod(teff,logg,microt,metal,outfile):
    path = os.path.abspath(os.path.dirname(__file__))
    kpath= path + '/modelatmospheres/'
    availteff=np.append(np.append(np.asarray([3500+(250*i) for i in range(35)]),np.asarray([12500+(500*i) for i in range(16)])),np.asarray([21000+(1000*i) for i in range(10)]))
    availlogg=np.asarray([0.5*i for i in range(11)])
    availmetal=np.asarray([-5,-4.5,-4,-3.5,-3,-2.8,-2.5,-2.3,-2.0,-1.8,-1.5,-1.3,-1.0,-0.8,-0.5,-0.3,0,0.3,0.5,0.8,1,1.5])
    atmotype='odfnew'
    ntau = 72

    v1=np.where(np.abs(availteff-teff) <= 0.1)
    v2=np.where(np.abs(availlogg-logg) <= 0.001)
    v3=np.where(np.abs(availmetal-metal) <= 0.001)

    if (min(availteff) <= teff <= max(availteff)) and (min(availlogg) <= logg <= max(availlogg)) and (min(availmetal) <= metal <= max(availmetal)):
        
        #If no interpolation needed, get correct model
        if len(v1[0]) > 0 and len(v2[0]) > 0 and len(v3[0]) > 0:
            model = rd_kmod(availteff[v1[0]][0],availlogg[v2[0]][0],availmetal[v3[0]][0])[0]
            
        else:
            tm1=max(availteff[np.where(availteff <= teff)])
            lm1=max(availlogg[np.where(availlogg <= logg)])
            mm1=max(availmetal[np.where(availmetal <= metal)])
            tp1=min(availteff[np.where(availteff >= teff)])
            lp1=min(availlogg[np.where(availlogg >= logg)])
            mp1=min(availmetal[np.where(availmetal >= metal)])
            
            ncols=10
            grid = np.zeros((2,2,2,ncols))
            
            if tp1 != tm1:
                mapteff = (teff-tm1)/(tp1-tm1)
            else:
                mapteff = 0.5
            if lp1 != lm1:
                maplogg = (logg-lm1)/(lp1-lm1)
            else:
                maplogg = 0.5
            if mp1 != mm1:
                mapmetal = (metal-mm1)/(mp1-mm1)
            else:
                mapmetal = 0.5
                
            for i in range(1,9):
                if i==1: model=rd_kmod(tm1,lm1,mm1)[0]
                if i==2: model=rd_kmod(tm1,lm1,mp1)[0]
                if i==3: model=rd_kmod(tm1,lp1,mm1)[0]
                if i==4: model=rd_kmod(tm1,lp1,mp1)[0]
                if i==5: model=rd_kmod(tp1,lm1,mm1)[0]
                if i==6: model=rd_kmod(tp1,lm1,mp1)[0]
                if i==7: model=rd_kmod(tp1,lp1,mm1)[0]
                if i==8: model=rd_kmod(tp1,lp1,mp1)[0]
                
                rhox = np.array([model[i][0] for i in range(ntau)])
                kappaross = np.array([model[i][4] for i in range(ntau)])
                tauross=np.zeros(ntau)
                tauross[0]=rhox[0]*kappaross[0]
                for ii in range(1,ntau):
                    tauross[ii]=simpson(kappaross[0:ii],rhox[0:ii])

                if i==1:
                    model1=model
                    tauross1=tauross
                elif i==2:
                    model2=model
                    tauross2=tauross
                elif i==3:
                    model3=model
                    tauross3=tauross
                elif i==4:
                    model4=model
                    tauross4=tauross
                elif i==5:
                    model5=model
                    tauross5=tauross
                elif i==6:
                    model6=model
                    tauross6=tauross
                elif i==7:
                    model7=model
                    tauross7=tauross
                elif i==8:
                    model8=model
                    tauross8=tauross
                    
            model = np.zeros((ntau,ncols))
            
            tauross = tauross1
            bot_tauross = min([tauross1[-1],tauross2[-1],tauross3[-1],tauross4[-1],tauross5[-1],tauross6[-1],tauross7[-1],tauross8[-1]])
            top_tauross = max([tauross1[0],tauross2[0],tauross3[0],tauross4[0],tauross5[0],tauross6[0],tauross7[0],tauross8[0]])
            
            tauross_new = np.interp(np.array(range(ntau)),np.array(range(ntau))[np.where((tauross >= top_tauross) & (tauross <= bot_tauross))], tauross[np.where((tauross >= top_tauross) & (tauross <= bot_tauross))])
            myinterp=interp1d(np.array(range(ntau))[np.where((tauross >= top_tauross) & (tauross <= bot_tauross))], tauross[np.where((tauross >= top_tauross) & (tauross <= bot_tauross))], fill_value='extrapolate',kind='cubic')
            tauross_new = myinterp(range(ntau))
                        
            for i in range(ntau):
                for j in range(ncols):
                    
                    interpmin = max(0,i-1000)
                    interpmax = min(i+1000,ntau-1)
                    
                    myinterp=interp1d(tauross1[interpmin:interpmax],model1[interpmin:interpmax,j], fill_value='extrapolate',kind='linear')
                    grid[0,0,0,j]=myinterp(tauross_new[i])
                    myinterp=interp1d(tauross2[interpmin:interpmax],model2[interpmin:interpmax,j], fill_value='extrapolate',kind='linear')
                    grid[0,0,1,j]=myinterp(tauross_new[i])
                    myinterp=interp1d(tauross3[interpmin:interpmax],model3[interpmin:interpmax,j], fill_value='extrapolate',kind='linear')
                    grid[0,1,0,j]=myinterp(tauross_new[i])
                    myinterp=interp1d(tauross4[interpmin:interpmax],model4[interpmin:interpmax,j], fill_value='extrapolate',kind='linear')
                    grid[0,1,1,j]=myinterp(tauross_new[i])
                    myinterp=interp1d(tauross5[interpmin:interpmax],model5[interpmin:interpmax,j], fill_value='extrapolate',kind='linear')
                    grid[1,0,0,j]=myinterp(tauross_new[i])
                    myinterp=interp1d(tauross6[interpmin:interpmax],model6[interpmin:interpmax,j], fill_value='extrapolate',kind='linear')
                    grid[1,0,1,j]=myinterp(tauross_new[i])
                    myinterp=interp1d(tauross7[interpmin:interpmax],model7[interpmin:interpmax,j], fill_value='extrapolate',kind='linear')
                    grid[1,1,0,j]=myinterp(tauross_new[i])
                    myinterp=interp1d(tauross8[interpmin:interpmax],model8[interpmin:interpmax,j], fill_value='extrapolate',kind='linear')
                    grid[1,1,1,j]=myinterp(tauross_new[i])
                                
                    myinterp2 = RegularGridInterpolator([[0,1],[0,1],[0,1]], grid[:,:,:,j])
                    model[i,j]=myinterp2([mapteff,maplogg,mapmetal])[0]
                                
        with open(outfile, 'w') as f:
            
            #Write the header
            teffstring="{:7.0f}".format(teff)                            
            loggstring="{:8.5f}".format(logg) 
            f.write('KURUCZ\n')
            f.write('TEFF' + teffstring + '.  GRAVITY' + loggstring + ' LTE \n')
            f.write('NTAU          72\n')
            
            # Write the body
            for i in range(ntau):
                f.write("{:15.8e}".format(model[i][0]).replace('e','E'))
                f.write("{:9.1f}".format(model[i][1]))
                f.write("{:10.3e}".format(model[i][2]).replace('e','E'))
                f.write("{:10.3e}".format(model[i][3]).replace('e','E'))
                f.write("{:10.3e}".format(model[i][4]).replace('e','E'))
                f.write("{:10.3e}".format(model[i][5]).replace('e','E'))
                f.write("{:10.3e}".format(model[i][6]).replace('e','E'))
                f.write("{:10.3e}".format(model[i][7]).replace('e','E'))
                f.write("{:10.3e}".format(model[i][8]).replace('e','E'))
                f.write("{:10.3e}".format(model[i][9]).replace('e','E'))
                f.write('\n')
                
            #Write the footer
            f.write(str(round(microt,2))+'\n')
            f.write('NATOMS         0	' + str(round(metal,1)) +'\n')

def py_kmod_main():
    
    path = os.path.abspath(os.path.dirname(__file__))
    if len(sys.argv) < 5:
        print("Usage: python pykmod.py (Teff (K)) (logg ([cm/s^2])) (vmicro (km/s)) ([Fe/H] (dex)) [output file path]")
        sys.exit()
    
    elif len(sys.argv) == 5:
        teff=float(sys.argv[1])
        logg=float(sys.argv[2])
        metal=float(sys.argv[4])
        microt=float(sys.argv[3])
        outfile=path + '/modelatmosphere.txt'
    
        pykmod(teff,logg,microt,metal,outfile)
        
    elif len(sys.argv) == 6:
        teff=float(sys.argv[1])
        logg=float(sys.argv[2])
        metal=float(sys.argv[4])
        microt=float(sys.argv[3])
        outfile=sys.argv[5]
        
        pykmod(teff,logg,microt,metal,outfile)

    elif len(sys.argv) > 6:
        print("Usage: python pykmod.py (Teff (K)) (logg ([cm/s^2])) (vmicro (km/s)) ([Fe/H] (dex)) [output file path]")
        sys.exit()

if __name__ == "__main__":
    py_kmod_main()            
