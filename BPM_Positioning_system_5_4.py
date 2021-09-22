# -*- coding: utf-8 -*-
"""
Created on Tue May  4 09:40:49 2021

@author: Xinru & Fangwei
"""

import MDAnalysis
from MDAnalysis import Universe
import numpy as np
import os

import logging
import sys
logging.basicConfig(
    stream=sys.stdout,
    level=logging.DEBUG,
    format='%(asctime)s %(name)s-%(levelname)s: %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S')


os.chdir('/pfs/proj/nobackup/fs/projnb10/snic2020-10-92/PlasmaMem/BPM2x2/python_scripts')

from class_BPM_Positioning_system_5_4 import Positioning_system

InputDir = '/pfs/proj/nobackup/fs/projnb10/snic2020-10-92/PlasmaMem/BPM2x2/download/'

OutputDir = '/pfs/proj/nobackup/fs/projnb10/snic2020-10-92/PlasmaMem/BPM2x2/download/'

mem = 1

U = Universe(InputDir+'1/conf.8.gro')

U2 = Universe(InputDir+'1/conf.8.gro')

start_run = 1

end_run = 60

Positioning_mem1 = np.zeros((end_run-start_run+1,4))

for run in range(start_run,end_run+1):
    
    try:
        pore = Positioning_system(mem,run,U,U2,InputDir)
        
        pore.coarse_search()  
        
        Time, Location = pore.search() 
        
        print(run,Time)
    
        print("The location of the pore:","(",Location[0],Location[1],")")
        
        Positioning_mem1[run-start_run] = np.array((run,Time,Location[0],Location[1]))
        
    except:
        
        print("run",run,"Error PASS")
        
        Positioning_mem1[run-start_run] = np.array((run,0,0,0))
    
        pass    

os.chdir(OutputDir+'1')

np.savetxt('positioning_mem1_5_4.txt',Positioning_mem1,delimiter=',')  

    
mem = 2

U = Universe(InputDir+'2/conf.8.gro')

U2 = Universe(InputDir+'2/conf.8.gro')

start_run = 1

end_run = 60

Positioning_mem2 = np.zeros((end_run-start_run+1,4))

for run in range(start_run,end_run+1):
    
    try:
    
        pore = Positioning_system(mem,run,U,U2,InputDir)
        
        pore.coarse_search()  
        
        Time, Location = pore.search() 
        
        print(run,Time)
    
        print("The location of the pore:","(",Location[0],Location[1],")")
        
        Positioning_mem2[run-start_run] = np.array((run,Time,Location[0],Location[1]))
        
    except:
        
        print("run",run,"Error PASS")
        
        Positioning_mem2[run-start_run] = np.array((run,0,0,0))
    
        pass    
    
os.chdir(OutputDir+'2')

np.savetxt('positioning_mem2_5_4.txt',Positioning_mem2,delimiter=',')  



mem = 3

U = Universe(InputDir+'3/conf.8.gro')

U2 = Universe(InputDir+'3/conf.8.gro')

start_run = 1

end_run = 60

Positioning_mem3 = np.zeros((end_run-start_run+1,4))

for run in range(start_run,end_run+1):
    
    try:
    
        pore = Positioning_system(mem,run,U,U2,InputDir)
        
        pore.coarse_search()  
        
        Time, Location = pore.search() 
        
        print(run,Time)
    
        print("The location of the pore:","(",Location[0],Location[1],")")
        
        Positioning_mem3[run-start_run] = np.array((run,Time,Location[0],Location[1]))
        
    except:
        
        print("run",run,"Error PASS")
        
        Positioning_mem3[run-start_run] = np.array((run,0,0,0))
    
        pass    
    
os.chdir(OutputDir+'3')

np.savetxt('positioning_mem3_5_4.txt',Positioning_mem3,delimiter=',')  



mem = 4

U = Universe(InputDir+'4/conf.8.gro')

U2 = Universe(InputDir+'4/conf.8.gro')

start_run = 1

end_run = 60

Positioning_mem4 = np.zeros((end_run-start_run+1,4))

for run in range(start_run,end_run+1):
    
    try:
        
        pore = Positioning_system(mem,run,U,U2,InputDir)
        
        pore.coarse_search()  
        
        Time, Location = pore.search() 
        
        print(run,Time)
    
        print("The location of the pore:","(",Location[0],Location[1],")")
        
        Positioning_mem4[run-start_run] = np.array((run,Time,Location[0],Location[1]))
        
    except:
        
        print("run",run,"Error PASS")
        
        Positioning_mem4[run-start_run] = np.array((run,0,0,0))
    
        pass    
    
os.chdir(OutputDir+'4')

np.savetxt('positioning_mem4_5_4.txt',Positioning_mem4,delimiter=',')  