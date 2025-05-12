# -*- coding: utf-8 -*-
"""
Created on Tue May  6 23:54:24 2025

@author: 17066
"""

import numpy as np
import shapely
from itertools import product
from copy import deepcopy
from collections import Counter
from smeshtools.hextools import *
from shapely import Polygon,GeometryCollection
from vectorScribe.objects.scribes import Template,Data_scribe
from tqdm import tqdm

class _Hex_mesh:
    """
    Base Object for working with hexagonal meshes.
    
    """
    
    def set_univ_vecs(self,**kwargs):
        self.univ_scribe.set_vecs(**kwargs)
    
    def set_univ_grps(self,**kwargs):
        self.univ_scribe.set_vecs(**kwargs)
    
    def set_mat_vecs (self,**kwargs):
        self.mat_scribe .set_vecs(**kwargs)
    
    def set_mat_grps (self,**kwargs):
        self.mat_scribe .set_grps(**kwargs)
    
    def set_elem_vecs(self,**kwargs):
        self.elem_scribe.set_vecs(**kwargs)
    
    def set_elem_grps(self,**kwargs):
        self.elem_scribe.set_grps(**kwargs)
    
    def constrain(self,origin,P,dz):
        """
        
        Parameters
        ----------
        origin : 1d iterable of float (length of 3)
            gives the centroid of the leftmost, lowermost hexagonal element
        P : float
            flat-to-flat pitch of the hexagonal array
        dz : 1d iterable of float
            gives the axial thickness of each aial layer in the hexagonal array
        
        """
        
        self.origin = origin
        self.P = P
        self.t = np.sqrt(3)/3*P
        self.dz = dz
        self.constrained = True
    
    def gen_interfaces(self):
        """
        
        Calculates all the interfaces on the hexagonal mesh
        
        """
        
        # Mesh Parameters
        nu = self.nu
        nv = self.nv
        nz = self.nz
        
        ele_void_id = self.ele_void_id
        ele_id_arr  = self.ele_id_arr
        
        # interfaces dictionary
        interfaces = []
        
        # Iterate through each Location
        us = np.linspace(0,nu-1,nu).astype(int)
        vs = np.linspace(0,nv-1,nv).astype(int)
        zs = np.linspace(0,nz-1,nz).astype(int)
        for z,v,u in tqdm(product(zs,vs,us)):
            ele_id = ele_id_arr[u,v,z]
            
            if ele_id==ele_void_id:
                continue
            else:
                interfaces = interface_calc(u,v,z,interfaces,self)
        
        self.interfaces = interfaces

class Hex_3D(_Hex_mesh):
    """
    Object for working with 3D hexagonal meshes
    
    """
    pass

class Hex_channel(_Hex_mesh):
    """
    Object for working with hexagonal meshes composed of "channels"
    
    """
    
    def __init__(self, channel_file, inp_id_map, channel_ids, out_id_map=None,
                 hex_type="x", channel_void_id="VD", inp_void_id="VD",
                 out_void_id="VD", ele_void_id="VD"):
        """
        
        Parameters
        ----------
        channel_file : str
            Name of the file that specifies the location of each channel in the
            xy plane
        
        inp_id_map : (m,n) iterable of str
            Provides universe composition in each channel. m channels each with
            n axial layers
        
        channel_ids : 1d iterable of str
            Provides the channel ids listed in ``channel-file`` and in order of
            the ``univ_id_map``
        
        mat_id_map : (m,n) iterable of str
            Provides material composition in each channel. m channels each with
            n axial layers. This is optional and intended for cases either when
            universe and material compositions differ or when universe and 
            material ids differ. 
        
        hex_type : str
            Specifies whether ``univ_file`` (and possibly ``mat_file``) 
            provides data in x-type or y-type hexagonal geometry
        
        channel_void_id : str
            Specifies the id for void regions in channel_file.
        
        univ_void_id : str
            Specifies the id for void regions in univ_id_map
        
        mat_void_id : str
            Specifies the id for void regions in mat_id_map
        
        """
        
        # Read in the hexagonal mesh
        channel_id_arr = hex_text_to_arr(channel_file,hex_type=hex_type,ndim=2)
        
        # Get dimensions
        nu,nv = np.shape(channel_id_arr)
        nc,nz = np.shape(inp_id_map)
        
        # if out_id_map is not specified, it is assumed to be the same as
        # inp_id_map
        if out_id_map is None:
            out_id_map  = deepcopy(inp_id_map)
            out_void_id = deepcopy(inp_void_id)
        else:
            # Need to check that univ_id_map and mat_id_map have the same 
            # dimensions
            
            # Need to check that univ_id_map and mat_id_map have voids in the same
            # locations
            
            pass
        
        # Obtain the id vectors
        inp_id_vec = np.unique(inp_id_map)
        out_id_vec = np.unique(out_id_map)
        
        # Build the id arrays
        inp_id_arr = np.full([nu,nv,nz],inp_void_id,dtype='<U7')
        out_id_arr = np.full([nu,nv,nz],out_void_id,dtype='<U7')
        for i,channel_id in enumerate(channel_ids):
            nds = np.where(channel_id_arr==channel_id)
            inp_id_arr[nds[0],nds[1],:] = inp_id_map[i,:]
            out_id_arr[nds[0],nds[1],:] = out_id_map[i,:]
        
        # Build the element id array
        vmask = inp_id_arr==inp_void_id
        umask = np.logical_not(vmask)
        ulocs = np.where(umask.transpose([2,1,0]))
        ne = np.count_nonzero(umask)
        
        ele_id_vec = np.linspace(0,ne-1,ne).astype(int).astype(str)
        ele_id_arr = np.full([nu,nv,nz],ele_void_id,dtype='<U7')
        ele_id_arr = ele_id_arr.transpose([2,1,0])
        ele_id_arr[ulocs] = ele_id_vec
        ele_id_arr = ele_id_arr.transpose([2,1,0])
        
        # Build the templates
        inp_template = Template(inp_id_vec, {"arr":inp_id_arr})
        out_template = Template(out_id_vec, {"arr":out_id_arr})
        ele_template = Template(ele_id_vec, {"arr":ele_id_arr})
        
        # Build the scribes
        inp_scribe = Data_scribe()
        inp_scribe.set_template(inp_template)
        out_scribe = Data_scribe()
        out_scribe.set_template(out_template)
        ele_scribe = Data_scribe()
        ele_scribe.set_template(ele_template)
        
        # Save data structures
        self.channel_void_id = channel_void_id
        self.inp_void_id     = inp_void_id
        self.out_void_id     = out_void_id
        self.ele_void_id     = ele_void_id
        self.channel_id_arr  = channel_id_arr
        self.inp_id_arr      = inp_id_arr
        self.out_id_arr      = out_id_arr
        self.ele_id_arr      = ele_id_arr
        
        self.nu = nu
        self.nv = nv
        self.nc = nc
        self.nz = nz
        
        self.inp_scribe  = inp_scribe
        self.out_scribe  = out_scribe
        self.ele_scribe  = ele_scribe
        self.hex_type    = hex_type
        self.constrained = False
    
    def gen_channel_centroids(self):
        """
        Generates the centroids in the x,y plane for each channel
        
        """
        
        # Get array size
        hex_type   = self.hex_type
        ch_id_arr  = self.channel_id_arr
        ch_void_id = self.channel_void_id
        
        nu     = self.nu
        nv     = self.nv
        origin = self.origin
        P      = self.P
        
        # Preallocate the channel centroid dictionary
        ch_centroids = {}
        for ch_id in np.unique(ch_id_arr):
            
            # Skip voids
            if ch_id==ch_void_id:
                continue
            
            ch_centroids[ch_id] = np.zeros([0,2])
        
        # Calculate spacing between hexagons
        if hex_type=="x":
            ux = P
            uy = 0
            vx = P/2
            vy = np.sqrt(3)*P/2
        
        elif hex_type=="y":
            ux = 0
            uy = P
            vx = np.sqrt(3)*P/2
            vy = P/2
        
        # Fill the channel centroid dictionary
        for v in range(nv):
            for u in range(nu):
                # Get channel id
                ch_id = ch_id_arr[u,v]
                
                # Skip voids
                if ch_id==ch_void_id:
                    continue
                
                # Hexagon centroid
                x0 = origin[0] + u*ux + v*vx
                y0 = origin[1] + u*uy + v*vy
                cen = [x0,y0]
                ch_centroids[ch_id] = np.vstack([ch_centroids[ch_id],cen])
        
        # Save values
        self.channel_centroids = ch_centroids
    
    def gen_channel_GeometryCollection(self):
        """
        Generates a Shapely GeometryCollection object to represent the 
        channel locations

        """
        
        # Get useful attributes
        channel_id_arr    = self.channel_id_arr
        channel_void_id   = self.channel_void_id
        channel_centroids = self.channel_centroids
        
        P = self.P
        t = self.t
        
        # Generate channel centroids
        self.gen_channel_centroids()
        
        # Create the channel_vec
        channel_id_vec = np.unique(channel_id_arr)
        mask           = channel_id_vec!=channel_void_id
        channel_id_vec = channel_id_vec[mask]
        
        # Point offsets from centroid of hexagon
        offset = np.array([[ 0  , t  ],
                           [ P/2, t/2],
                           [ P/2,-t/2],
                           [ 0  ,-t  ],
                           [-P/2,-t/2],
                           [-P/2, t/2]])
        
        # Create the unioned geometries
        geom_vec = []
        for channel_id in channel_id_vec:
            # Grab the centroids corresponding to a given channel
            centroids = channel_centroids[channel_id]
            hexes = []
            
            # Create a hexagon for each centroid in the list
            for centroid in centroids:
                x0,y0 = centroid
                hex_ = deepcopy(offset)
                hex_[:,0] += x0
                hex_[:,1] += y0
                
                hex_ = Polygon(hex_)
                hexes.append(hex_)
            
            # Create a geometry collection from all the hexes
            geom = shapely.union_all(hexes)
            geom_vec.append(geom)
        
        # Save the geometry collection
        self.channel_GeometryCollection = GeometryCollection(geom_vec)

### HELPER FUNCTIONS AND OBJECTS ##############################################

class Interface:
    """Useful object for storing data about a given interface"""
    
    def __init__(self, inp_ids, out_ids, ele_ids, A, dlts):
        """
        
        Parameters
        ----------
        inp_ids : iterable of str of length 2
            provides the input ids on both sides of the interface
        out_ids : iterable of str of length 2
            provides the output ids on both sides of the interface
        ele_ids : iterable of str of length 2
            provides the element ids on both sides of the interface
        A : float
            area of the interface
        dlts : iterable of float of length 2
            provides the distance between the centroid and the interface on 
            both sides of the interface
        
        """
        
        # need to check that the length of each of these attribures is 2
        
        self.inp_ids = inp_ids
        self.out_ids = out_ids
        self.ele_ids = ele_ids
        self.A = A
        self.dlts = dlts
    
    def __eq__(self, o):
        """
        
        Two interfaces are equivalent if both contain the same element ids.
        
        """
        return Counter(self.ele_ids)==Counter(o.ele_ids)

def append_interface(interface_bool, interface_type, ndx_a, ndx_b, interfaces,
                     obj):
    """
    Function that creates and appends a given interface to the interfaces list
    
    Parameters
    ----------
    interface_bool : boolean
        Indicates wheter element a interfaces with another element
    interface_type : str
        Either "base" or "side" indicates whether the interface is on a base
        or side of element a
    ndx_a : 1d iterable of int of length 3
        Provides the mesh index of element a
    ndx_b : 1d iterable of int of length 3
        Provides the mesh index of element b
    interfaces : iterable of Interface objects
        The list of interfaces to append new interfaces to
    obj : Hex_channel object
        The Hex_channel object from which to pull relevant attributes
    
    """
    
    # Extract useful attributes
    inp_void_id = obj.inp_void_id
    out_void_id = obj.out_void_id
    ele_void_id = obj.ele_void_id
    
    inp_id_arr = obj.inp_id_arr
    out_id_arr = obj.out_id_arr
    ele_id_arr = obj.ele_id_arr
    
    P  = obj.P
    t  = obj.t
    dz = obj.dz
    
    # Calculate side a parameters
    inp_id_a = inp_id_arr[ndx_a]
    out_id_a = out_id_arr[ndx_a]
    ele_id_a = ele_id_arr[ndx_a]
    
    if interface_type=="base":
        A = 3*np.sqrt(3)/2*t**2
        dlt_a = dz[ndx_a[2]]
    elif interface_type=="side":
        A = dz[ndx_a[2]]*t
        dlt_a = P/2
    
    # Get id on side b
    if interface_bool:
        # Side b id is another element
        inp_id_b = inp_id_arr[ndx_b]
        out_id_b = out_id_arr[ndx_b]
        ele_id_b = ele_id_arr[ndx_b]
    
    else:
        # Side b id is a void
        inp_id_b = inp_void_id
        out_id_b = out_void_id
        ele_id_b = ele_void_id
    
    # Set the ids for the interface
    inp_ids = [inp_id_a,inp_id_b]
    out_ids = [out_id_a,out_id_b]
    ele_ids = [ele_id_a,ele_id_b]
    
    # Update interface boolean
    if ele_id_b==ele_void_id:
        interface_bool=False
    
    # Create the interface
    if interface_bool:
        
        if interface_type=="base":
            dlt_b = dz[ndx_b[2]]
        elif interface_type=="side":
            dlt_b = P/2
        
        dlts = [dlt_a,dlt_b]
        
        interface = Interface(inp_ids,out_ids,ele_ids,A,dlts)
        
        if not np.isin(interface,interfaces):
            interfaces.append(interface)
        else:
            pass
    
    else:
        
        dlts = [dlt_a,2/3]
        
        interface = Interface(inp_ids,out_ids,ele_ids,A,dlts)
        interfaces.append(interface)
    
    return interfaces

def interface_calc(u, v, z, interfaces, obj):
    """
    Creates and appends all interfaces associated with element "a", the 
    element located at the mesh index (u,v,z)
    
    Parameters
    ----------
    u : int
        mesh index along the u direction
    v : int
        mesh index along the v direction
    z : int
        mesh index along the z direction
    interfaces : iterable of Interface objects
        The list of interfaces to append new interfaces to
    obj : Hex_channel object
        The Hex_channel object from which to pull relevant attributes
    
    Returns
    -------
    interfaces : iterable of Interface objects
        The list of interfaces with new interfaces appended

    """
    
    # Extract useful attributes
    nu = obj.nu
    nv = obj.nv
    nz = obj.nz
    
    # Extract data for side a
    ndx_a = (u,v,z)
    
    # Preallocate interface booleans
    a_bool = True
    b_bool = True
    c_bool = True
    d_bool = True
    e_bool = True
    f_bool = True
    g_bool = True
    h_bool = True
    
    # Check which faces on element "a" interface with other elements
    if u==0:
        d_bool = False
        f_bool = False
    if u==nu-1:
        c_bool = False
        e_bool = False
    if v==0:
        b_bool = False
        c_bool = False
    if v==nv-1:
        f_bool = False
        g_bool = False
    if z==0:
        a_bool = False
    if z==nz-1:
        h_bool = False
    
    if u==22:
        pass
    
    # Append all interfaces
    interfaces = append_interface(a_bool,"base",ndx_a,(u  ,v  ,z-1),
                                  interfaces,obj)
    interfaces = append_interface(b_bool,"side",ndx_a,(u  ,v-1,z  ),
                                  interfaces,obj)
    interfaces = append_interface(c_bool,"side",ndx_a,(u+1,v-1,z  ),
                                  interfaces,obj)
    interfaces = append_interface(d_bool,"side",ndx_a,(u-1,v  ,z  ),
                                  interfaces,obj)
    interfaces = append_interface(e_bool,"side",ndx_a,(u+1,v  ,z  ),
                                  interfaces,obj)
    interfaces = append_interface(f_bool,"side",ndx_a,(u-1,v+1,z  ),
                                  interfaces,obj)
    interfaces = append_interface(g_bool,"side",ndx_a,(u  ,v+1,z  ),
                                  interfaces,obj)
    interfaces = append_interface(h_bool,"base",ndx_a,(u  ,v  ,z+1),
                                  interfaces,obj)
    
    return interfaces