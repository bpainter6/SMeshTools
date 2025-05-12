# -*- coding: utf-8 -*-
"""
Created on Tue May  6 23:51:32 2025

@author: 17066
"""

import numpy as np
from copy import deepcopy
from hextools import hex_text_to_arr
from shapely import Polygon,GeometryCollection

class Hex_3D():
    
    def __init__(self,filename,hex_type="x"):
        """
        Parameters
        ----------
        filename : str
            Name of the file that provides the geometry of the hexagonal
            mesh

        """
        
        # Read in the hexagonal mesh
        self.id_arr = hex_text_to_arr(filename,hex_type=hex_type,ndim=3)
        self.hex_type = hex_type
        self.constrained = False
        
        # Create the channel map
        nu,nv,nz = np.shape(self.arr)
        nc = nu*nv # Number of channels
        self.nu = nu
        self.nv = nv
        self.nz = nz
        self.nc = nc
        self.map = self.arr.transpose([1,0,2]).reshape(nu*nv,nz)