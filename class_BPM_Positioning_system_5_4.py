# -*- coding: utf-8 -*-
"""
Created on Tue May  4 09:39:30 2021

@author: Xinru & Fangwei
"""

import numpy as np
import MDAnalysis
from MDAnalysis import Universe
import numpy as np
import os

class Positioning_system:
        
    
    def  __init__(self,mem,run,U,U2,InputDir):
        
        self.mem = mem
        
        self.run = run
        
        self.dir = InputDir
        
        self.U = U
        
        self.U2 = U2
        
        self.coarse_frame = 6   
        
        self.region_i = 0
        
        self.region_j = 0
        
        self.nodes = [6*i for i in range(2,11)]
        
        self.coarse_search_stop = 70
        
        self.search_after_60 = False
        
        self.time = 0
        
        self.location = np.array([0,0,0])
        
        
    def _check_next_frame(self,frame,PW):
        
         self.U2 =  self.U2.load_new(self.dir+str(self.mem)+'/'+str(self.mem)+'-'+str(self.run)+'/conf.9.xtc')
        
         x = str(PW.positions[0][0])
                
         y = str(PW.positions[0][1])
                
         z = str(PW.positions[0][2])
        
         for TS in self.U2.trajectory:
            
            if (TS.frame == frame+1):
            
                NF = self.U2.select_atoms("(name W and point "+x+" "+y+" "+z+" 15)")
            
            if (TS.frame == frame+2):
            
                NNF = self.U2.select_atoms("(name W and point "+x+" "+y+" "+z+" 15)")
                
                return NF,NNF
           

            
    def _check_around(self,frame,PW):
        
        self.U2 =  self.U2.load_new(self.dir+str(self.mem)+'/'+str(self.mem)+'-'+str(self.run)+'/conf.9.xtc')
        
        for Ts in self.U2.trajectory:
            
            if Ts.frame == frame:
                
                Up = self.U2.select_atoms("(name W and prop z >= "+str(PW.positions[1][2])+" and prop z <= "+str(PW.positions[1][2]+12)+" and prop x >= "+str(PW.positions[1][0]-12)+" and prop x <= "+str(PW.positions[1][1]+12)+" and prop y >= "+str(PW.positions[1][1]-12)+" and prop y <= "+str(PW.positions[1][1]+12)+")",updating = True)
    
                Down = self.U2.select_atoms("(name W and prop z >= "+str(PW.positions[1][2]-12)+" and prop z <= "+str(PW.positions[1][2])+" and prop x >= "+str(PW.positions[1][0]-12)+" and prop x <= "+str(PW.positions[1][1]+12)+" and prop y >= "+str(PW.positions[1][1]-12)+" and prop y <= "+str(PW.positions[1][1]+12)+")",updating = True)
    
                return Up, Down
        


    def coarse_search(self):
                
        self.U =  self.U.load_new(self.dir+str(self.mem)+'/'+str(self.mem)+'-'+str(self.run)+'/conf.9.xtc')
        
        for ts in  self.U.trajectory:
                       
            if ts.frame > self.coarse_search_stop:
                
                self.search_after_60 = True
                
                return self
            
            if ts.frame < self.coarse_frame:
                
                continue
                        
            for node in self.nodes:
                
                if ts.frame != node:
                    
                    continue
                
                elif ts.frame == node:
                    
                    for i in range (0,10):
                
                        x1 = str(i*30)
                        
                        if i == 9:
                                
                            x2 = str(300)
                                
                        else:
                        
                            x2 = str(i*30+40)
                            
                        for j in range (0,10):
                                
                            y1 = str(j*30)
                            
                            if j == 9:
                                
                                y2 = str(300)
                                
                            else:
                        
                                y2 = str(j*30+40)
                                
                            try:
                            
                                u =  self.U.select_atoms("(all and prop x>= "+x1+" and prop x<="+x2+" and prop y>= "+y1+" and prop y<="+y2+")", updating=True)
                                
                                com = (u.select_atoms("(all and (not resname PW ION))", updating=True)).center_of_mass()                           
                              
                                if ts.frame <= self.coarse_search_stop:
                                    
                                    com_upper = (u.select_atoms("(name PO4 and prop z > "+str(com[2])+")", updating=True)).center_of_mass()
                                
                                    com_lower = (u.select_atoms("(name PO4 and prop z < "+str(com[2])+")", updating=True)).center_of_mass()
                            
                                    thickness = com_upper[2] - com_lower[2]
                                
                                    dz = thickness/5
                            
                                    PW = u.select_atoms("(name W and prop z >= "+str(com[2]-dz/2)+" and prop z <= "+str(com[2]+dz/2)+")", updating=True)
                                     
                                    if len(PW) > 12:
                                        
                                        self.coarse_frame = ts.frame
                                        
                                        self.region_i = i
                                        
                                        self.region_j = j 
                        
                                        return self
                                    
                            except:
                                        
                                print(ts.frame)
                                print('Error PASS!!!')
                                        
                                pass
         
    
    def search(self):
        
        if self.search_after_60:
            
            self._Search_After()
            
        else:
            
            self._Precise_search()
            
            if self.search_after_60:
                
                self._Search_After()
            
        return self.time, self.location

        
        
        
       
    def _Search_After(self):
    
         self.U =  self.U.load_new(self.dir+str(self.mem)+'/'+str(self.mem)+'-'+str(self.run)+'/conf.9.xtc')
        
         print("search_after")         
         
         for ts in self.U.trajectory:
                
            if ts.frame <= self.coarse_search_stop:
                    
                continue
            
            print(ts.frame)
                
            for i in range (0,10):
                    
                x1 = str(i*30)
                        
                if i == 9:
                                
                    x2 = str(300)
                                
                else:
                        
                    x2 = str(i*30+40)
                            
                for j in range (0,10):
                                    
                    y1 = str(j*30)
                                
                    if j == 9:
                                    
                        y2 = str(300)
                                    
                    else:
                            
                        y2 = str(j*30+40)
                                
                    u =  self.U.select_atoms("(all and prop x>= "+x1+" and prop x<="+x2+" and prop y>= "+y1+" and prop y<="+y2+")", updating=True)
                                    
                    com = (u.select_atoms("(all and not resname PW ION)", updating=True)).center_of_mass()
                    
                    try:
                         
                        if ts.frame >= 60 and ts.frame <= 70:
                                        
                            com_upper = (u.select_atoms("(name PO4 and prop z > "+str(com[2])+")", updating=True)).center_of_mass()
                                        
                            com_lower = (u.select_atoms("(name PO4 and prop z < "+str(com[2])+")", updating=True)).center_of_mass()
                                    
                            thickness = com_upper[2] - com_lower[2]
                                        
                            dz = thickness/5
                                    
                            PW = u.select_atoms("(name W and prop z >= "+str(com[2]-dz/2)+" and prop z <= "+str(com[2]+dz/2)+")", updating=True)
                                                                                    
                            if len(PW) > 6:
                                                                                    
                                NF, NNF = self._check_next_frame(ts.frame,PW)
                                            
                                if len(NF)>0 and len(NNF)>0:
                                                
                                    Up, Down = self._check_around(ts.frame,PW)
                                                
                                    if len(Up)>0 and len(Down)>0: 
                                                    
                                        pore_location = PW.center_of_mass() 
                                        
                                        self.time = ts.frame
                                        
                                        self.location = pore_location
                                        
                                        return self
                                
                        if ts.frame > 70:
                                    
                             print('NO PORE BEFORE 70!!!')
                                    
                             return self
        
                    except:
                                        
                        print(ts.frame)
                        print('Error PASS!!!')
                                        
                        pass
             
         
            
    def _Precise_search(self):
        
        print("precise_search")
        
        self.U =  self.U.load_new(self.dir+str(self.mem)+'/'+str(self.mem)+'-'+str(self.run)+'/conf.9.xtc')
        
        for ts in self.U.trajectory:
                                 
            if ts.frame <= self.coarse_frame-6:
                    
                continue
            
            print(ts.frame)
            
            if ts.frame > self.coarse_frame:
                
                self.search_after_60 = True
                
                self.coarse_search_stop = self.coarse_frame
                
                return self
                    
            for i in range (self.region_i-2,self.region_i+2):
                
                if i < 0:
                        
                    i = 10+i
                    
                if i > 9:
                    
                    i = i-10
                    
                x1 = str(i*30)
                        
                if i == 9:
                            
                    x2 = str(300)
                            
                else:
                    
                    x2 = str(i*30+40)
                            
                for j in range (self.region_j-2,self.region_j+2):
                            
                    if j < 0:
                        
                        j = 10+j
                        
                    if j > 9:
                    
                        j = j-10
                        
                    y1 = str(j*30)
                        
                    if j == 9:
                            
                        y2 = str(300)
                            
                    else:
                    
                        y2 = str(j*30+40)
                        
                    u = self.U.select_atoms("(all and prop x>= "+x1+" and prop x<="+x2+" and prop y>= "+y1+" and prop y<="+y2+")", updating=True)
                    
                    com = (u.select_atoms("(all and not resname PW ION)", updating=True)).center_of_mass()
                    
                    PW = u.select_atoms("(name W and prop z >= "+str(com[2]-5)+" and prop z <= "+str(com[2]+5)+")", updating=True)
                    
                           
                    if len(PW) > 6:
                                                                    
                        NF, NNF = self._check_next_frame(ts.frame,PW)
                        
                        if len(NF)>0 and len(NNF)>0:
                                        
                            #Up, Down = self._check_around(ts.frame,PW)
                                        
                            #if len(Up)>0 and len(Down)>0:
                                        
                            #pore_location = PW.center_of_mass()
                                    
                            #return i,j,ts.frame,PW,NF,Up
                                        
                            pore_location = PW.center_of_mass()
                            
                            self.time = ts.frame
                                    
                            self.location = pore_location
                                            
                            return self