# -*- coding: utf-8 -*-
"""
Created on Wed Apr  2 22:23:57 2025

@author: 17066
"""

import numpy as np
import shapely
from shapely import Polygon,GeometryCollection
from copy import deepcopy

def rect_gc(xmin,xmax,ymin,ymax,channel_map,void_id="VD"):
    """Creates a GeometryCollection representing a rectangular mesh
    
    Parameters
    ----------
    channel_map - 2d iterable of str
        provides the channel id for each element of the mesh. Channels 
        appearing multiple times are unioned into a MultiPolygon object
    void_id - str
        indicates what regions in the mesh are skipped when creating a
        (multi)polygon
    
    """
    
    # Get dimensions of the mesh
    ny,nx = np.shape(channel_map)
    
    # Create the channel_vec
    channel_vec = np.unique(channel_map)
    mask = channel_vec!=void_id
    channel_vec = channel_vec[mask]
    
    # Width of each element
    dx = (xmax-xmin)/nx
    dy = (ymax-ymin)/ny
    
    # Point offsets from centroid of rectangle
    offset = np.array([[ dx/2, dy/2],
                       [ dx/2,-dy/2],
                       [-dx/2,-dy/2],
                       [-dx/2, dy/2]])
    
    rects = np.zeros((ny,nx),dtype=object)
    for j in range(ny):
        for i in range(nx):
            x0 = xmin+i*dx+dx/2
            y0 = ymin+j*dy+dy/2
            
            rect = deepcopy(offset)
            rect[:,0] += x0
            rect[:,1] += y0
            
            rect = Polygon(rect)
            rects[j,i] = rect
    
    # Create the unioned geometries
    geom_vec = []
    for channel_id in channel_vec:
        mask = channel_map==channel_id
        geom = shapely.union_all(rects[mask])
        geom_vec.append(geom)
    
    return channel_vec,GeometryCollection(geom_vec)

def x_hex_gc(origin,P,channel_map,void_id="VD"):
    """Creates a GeometryCollection representing a hexagonal mesh"""
    
    # Get dimensions of the mesh
    ny,nx = np.shape(channel_map)
    
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
    hexes = np.zeros([ny,nx],dtype=object)
    for u in range(ny):
        for v in range(nx):
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

def y_hex_gc(origin,P,channel_map,void_id="VD"):
    """Creates a GeometryCollection representing a hexagonal mesh"""
    
    # Get dimensions of the mesh
    ny,nx = np.shape(channel_map)
    
    # Create the channel_vec
    channel_vec = np.unique(channel_map)
    mask = channel_vec!=void_id
    channel_vec = channel_vec[mask]
    
    # Side length
    t = np.sqrt(3)/3*P
    
    # Spacing between hexagons
    ux = np.sqrt(3)*P/2
    uy = P/2
    vx = 0
    vy = P
    
    # Point offsets from centroid of hexagon
    offset = np.array([[ t/2, P/2],
                       [ t  , 0  ],
                       [ t/2,-P/2],
                       [-t/2,-P/2],
                       [-t  , 0  ],
                       [-t/2, P/2]])
    
    # Create each hexagon
    hexes = np.zeros((ny,nx),dtype=object)
    for u in range(nx):
        for v in range(ny):
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