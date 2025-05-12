# -*- coding: utf-8 -*-
"""
Created on Mon May  5 23:15:42 2025

@author: 17066
"""

import numpy as np
from copy import deepcopy
from hextools import hex_text_to_arr
from shapely import Polygon,GeometryCollection

class Hex_2D():
    
    def __init__(self,filename,hex_type="x"):
        """
        Parameters
        ----------
        filename : str
            Name of the file that provides the geometry of the hexagonal
            mesh

        """
        
        # Read in the hexagonal mesh
        self.id_arr = hex_text_to_arr(filename,hex_type=hex_type,ndim=2)
        self.hex_type = hex_type
        self.constrained = False
        
        # Create the channel map
        nu,nv = np.shape(self.arr)
        self.nu = nu
        self.nv = nv
        self.nc = nc
        self.map = self.arr.transpose([1,0]).reshape(nu*nv)
    
    def constrain(self,origin,pitch):
        """

        Parameters
        ----------
        origin : 1d iterable of float
            centroid of the bottom left hexagon
        pitch : float
            flat-to-flat distance of the hexagonal array

        """
        
        self.origin = origin
        self.pitch = pitch
        self.constrained = True
    
    def gen_channel_centroids(self):
        """
        Generates the centroids in the x,y plane for each channel
        
        """
        
        # Get array size
        nu = self.nu
        nv = self.nv
        nc = self.nc
        
        # Get geometry of the mesh
        origin = self.origin
        P = self.pitch
        
        # Side length
        t = np.sqrt(3)/3*P
        
        if self.hex_type=="x":
            # Spacing between hexagons
            ux = P/2
            uy = np.sqrt(3)*P/2
            vx = P
            vy = 0
            
            # Create each hexagon
            centroids = np.zeros([nc,2],dtype=object)
            channel_ndx = 0
            for v in range(nv):
                for u in range(nu):
                    # Hexagon centroid
                    x0 = origin[0] + u*ux + v*vx
                    y0 = origin[1] + u*uy + v*vy
                    centroids[c,:] = np.array([x0,y0])
                    c += 1
        else:
            # Spacing between hexagons
            ux = np.sqrt(3)*P/2
            uy = P/2
            vx = 0
            vy = P
            
            # Create each hexagon
            centroids = np.zeros([nu*nv,2],dtype=object)
            c = 0
            for v in range(nv):
                for u in range(nu):
                    # Hexagon centroid
                    x0 = origin[0] + u*ux + v*vx
                    y0 = origin[1] + u*uy + v*vy
                    centroids[c,:] = np.array([x0,y0])
                    c += 1
        
        self.channel_centroids = centroids
    
    def gen_centroids(self):
        """
        Generates the centroids associated with each hexagon in the 
        hexagonal array
        
        """
        
        # Generate channel centroids
        self.gen_channel_centroids()
        
        # Get array size
        nu = self.nu
        nv = self.nv
        nc = self.nc
        
        # Preallocate centroids dictionary
        centroids = {}
        for Id in np.unique(self.arr):
            centroids[Id] = np.zeros([0,2])
        
        # Append centroids
        for c in nc:
            Id = self.map[c]
            centroid = self.channel_centroids[c,:].reshape(1,-1)
            centroids[Id] = np.vstack([centroids[Id],centroid])
    
    def gen_interfaces(self):
        nu,nv,nz = 
    
    def gen_GeometryCollection(self):
        """
        Generates a Shapely GeometryCoolection object to represent the 
        hexagonal array in the xy plane

        """
        # Get dimensions of the mesh
        nu,nv = np.shape(self.arr)
        
        # Get geometry of the mesh
        origin = self.origin
        P = self.pitch
        
        # Create the channel_vec
        channel_vec = np.unique(channel_map)
        mask = channel_vec!=void_id
        channel_vec = channel_vec[mask]
        
        # Side length
        t = np.sqrt(3)/3*P
        
        # Spacing between hexagons
        ux = P/2
        uy = np.sqrt(3)*P/2
        vx = P
        vy = 0
        
        # Point offsets from centroid of hexagon
        offset = np.array([[ 0  , t  ],
                           [ P/2, t/2],
                           [ P/2,-t/2],
                           [ 0  ,-t  ],
                           [-P/2,-t/2],
                           [-P/2, t/2]])
        
        # Create each hexagon
        hexes = np.zeros([nu,nv],dtype=object)
        for u in range(nu):
            for v in range(nv):
                # Hexagon centroid
                x0 = origin[0] + u*ux + v*vx
                y0 = origin[1] + u*uy + v*vy
                
                hex_ = deepcopy(offset)
                hex_[:,0] += x0
                hex_[:,1] += y0
                
                hex_ = Polygon(hex_)
                hexes[v,u] = hex_
        
        # Create the unioned geometries
        geom_vec = []
        for channel_id in channel_vec:
            mask = channel_map==channel_id
            geom = shapely.union_all(hexes[mask])
            geom_vec.append(geom)
        
        return channel_vec,GeometryCollection(geom_vec)
    
    def id_mat(self):
        pass
    
    