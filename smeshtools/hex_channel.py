# -*- coding: utf-8 -*-
"""
Created on Tue May  6 23:54:24 2025

@author: 17066
"""

import numpy as np
import shapely
from copy import deepcopy
from smeshtools.hextools import hex_text_to_arr
from shapely import Polygon,GeometryCollection
from vectorScribe.objects.scribes import Template,Data_scribe
from scipy.spatial import KDTree

class _Hex_3D:
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
        
        if not self.constrained:
            # Needs to raise error if not constrained
            pass
        
        # Extract useful attributes
        void_id = self.void_id
        
        P = self.P
        t = self.t
        
        # Get interface array centroids
        self.gen_interface_attributes()
        
        int_inp_id_flt = self.int_inp_id_flt
        int_out_id_flt = self.int_out_id_flt
        
        int_dz      = self.int_dz
        int_w_flt   = self.int_uvw_flt[:,2]
        int_xyz_flt = self.int_xyz_flt
        
        # Use the centroids to create a KDTree
        tree = KDTree(int_xyz_flt)
        
        # Query the KDTree for interfaces
        int_pairs = tree.query_pairs(r=1,eps=0.1)
        int_pairs = np.array(list(int_pairs))
        
        int_inp_id_pairs = flt_to_pairs(int_inp_id_flt, int_pairs, '<U7')
        
        # Filter out double void regions
        vmask0 = int_inp_id_pairs[:,0]==void_id
        vmask1 = int_inp_id_pairs[:,1]==void_id
        vmask  = np.logical_and(vmask0,vmask1)
        imask  = np.logical_not(vmask)
        
        int_pairs        = int_pairs[imask,:]
        int_inp_id_pairs = int_inp_id_pairs[imask,:]
        
        # Convert the remaining data
        int_out_id_pairs = flt_to_pairs(int_out_id_flt, int_pairs, '<U7')
        
        # Hexagonal index along double direction for each element
        int_w_pairs = flt_to_pairs(int_w_flt, int_pairs, int)
        
        smask = int_w_pairs[:,0]==int_w_pairs[:,1] # Side interface pairs
        bmask = np.logical_not(smask)              # Base interface pairs
        vmask = int_inp_id_pairs==void_id          # Void regions
        
        # Distances between centroids
        int_dlt_pairs = np.zeros_like(int_w_pairs,dtype=float)
        
        # Side interfaces
        int_dlt_pairs[smask,:] = P/2
        
        # Base interfaces
        int_dlt_pairs[bmask,:]  = flt_to_pairs(int_dz,int_w_pairs[bmask,:],
                                               float)
        int_dlt_pairs[bmask,:] /= 2
        
        # Void regions
        int_dlt_pairs[vmask] = 2/3
        
        # Interface area
        int_area_pairs = np.zeros_like(smask, dtype=float)
        
        # Side interfaces
        int_area_pairs[smask]  = int_dz[int_w_pairs[smask,0]]
        int_area_pairs[smask] *= t
        
        # Base interfaces
        int_area_pairs[bmask] = 3*np.sqrt(3)*t**2/2
        
        # Save results
        self.int_inp_id_pairs = int_inp_id_pairs
        self.int_out_id_pairs = int_out_id_pairs
        
        self.int_dlt_pairs  = int_dlt_pairs
        self.int_area_pairs = int_area_pairs

class Hex_mesh_3D(_Hex_3D):
    """
    Object for working with 3D hexagonal meshes
    
    """
    pass

class Hex_channel(_Hex_3D):
    """
    Object for working with hexagonal meshes composed of "channels"
    
    """
    
    def __init__(self, channel_file, inp_id_map, channel_ids, out_id_map=None,
                 hex_type="x", channel_void_id="VD", void_id="VD"):
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
        
        void_id : str
            Specifies the id for void regions
        
        """
        
        # Read in the hexagonal mesh
        channel_id_arr = hex_text_to_arr(channel_file,hex_type=hex_type,ndim=2)
        
        # Get dimensions
        nu,nv = np.shape(channel_id_arr)
        nc,nw = np.shape(inp_id_map)
        
        # if out_id_map is not specified, it is assumed to be the same as
        # inp_id_map
        if out_id_map is None:
            out_id_map  = deepcopy(inp_id_map)
        else:
            # Need to check that univ_id_map and mat_id_map have the same 
            # dimensions
            
            # Need to check that univ_id_map and mat_id_map have voids in the same
            # locations
            
            pass
        
        # Obtain the input/output id vectors
        inp_id_vec = np.unique(inp_id_map)
        out_id_vec = np.unique(out_id_map)
        
        # Build the input/output id flat and array
        inp_id_arr = np.full([nu,nv,nw], void_id, dtype='<U7')
        out_id_arr = np.full([nu,nv,nw], void_id, dtype='<U7')
        for i,channel_id in enumerate(channel_ids):
            nds = np.where(channel_id_arr==channel_id)
            inp_id_arr[nds[0],nds[1],:] = inp_id_map[i,:]
            out_id_arr[nds[0],nds[1],:] = out_id_map[i,:]
        
        # Build the templates
        inp_template = Template(inp_id_vec, {"arr":inp_id_arr})
        out_template = Template(out_id_vec, {"arr":out_id_arr})
        
        # Build the scribes
        inp_scribe = Data_scribe()
        inp_scribe.set_template(inp_template)
        out_scribe = Data_scribe()
        out_scribe.set_template(out_template)
        
        # Save data structures
        self.hex_type        = hex_type
        self.channel_void_id = channel_void_id
        self.void_id         = void_id
        self.channel_id_arr  = channel_id_arr
        
        self.nu = nu
        self.nv = nv
        self.nc = nc
        self.nw = nw
        
        self.inp_id_arr = inp_id_arr
        self.out_id_arr = out_id_arr
        
        self.inp_scribe = inp_scribe
        self.out_scribe = out_scribe
        
        self.constrained = False
    
    def gen_interface_attributes(self):
        """Generates centroid values for interface calculations"""
        
        if not self.constrained:
            # Needs to raise error if not constrained
            pass
        
        # Get useful attributes
        hex_type = self.hex_type
        void_id  = self.void_id
        nu       = self.nu
        nv       = self.nv
        nw       = self.nw
        
        inp_id_arr = self.inp_id_arr
        out_id_arr = self.out_id_arr
        
        dz = self.dz
        
        # Build the interface input/output array and flatten
        int_inp_id_arr = np.full([nu+2,nv+2,nw+2], void_id, dtype='<U7')
        int_out_id_arr = np.full([nu+2,nv+2,nw+2], void_id, dtype='<U7')
        int_inp_id_arr[ 1:-1, 1:-1, 1:-1] = inp_id_arr
        int_out_id_arr[ 1:-1, 1:-1, 1:-1] = out_id_arr
        int_inp_id_flt = int_inp_id_arr.transpose([2,1,0]).flatten()
        int_out_id_flt = int_out_id_arr.transpose([2,1,0]).flatten()
        
        # Uniformly spaced mesh for identifying pairs
        uni_dz  = np.ones(nw+2)
        int_res = centroid_calc((nu+2,nv+2,nw+2), (0,0,0), 1, uni_dz, hex_type)
        
        # Axial heights of the interface mesh
        int_dz = np.hstack([ dz[0], dz, dz[-1] ])
        
        # Save results
        self.int_inp_id_flt = int_inp_id_flt
        self.int_out_id_flt = int_out_id_flt
        
        self.int_dz = int_dz
        
        self.int_uvw_flt = int_res[2]
        self.int_xyz_flt = int_res[3]
    
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
        
        # Hexagonal indices
        u = np.linspace(0,nu-1,nu,dtype=int).reshape(-1, 1)
        v = np.linspace(0,nv-1,nv,dtype=int).reshape( 1,-1)
        
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
        
        # Centroids in cartesian coordinates as an array
        xa = origin[0] + ux*u + vx*v
        ya = origin[1] + uy*u + vy*v
        
        # Populatae the channel centroid dictionary
        ch_centroids = {}
        for ch_id in np.unique(ch_id_arr):
            
            # Skip voids
            if ch_id==ch_void_id:
                continue
            
            # Identify all locations of ch_id
            mask = ch_id_arr==ch_id
            
            # Grab centroids at those locations
            x0 = xa[mask].reshape(-1,1)
            y0 = ya[mask].reshape(-1,1)
            
            # Save the results
            ch_centroids[ch_id] = np.hstack([x0,y0])
        
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

def flt_to_pairs(dat_flt,pairs,typ):
    # Convert pair indices to inp/out ids
    dat_pairs = np.zeros_like(pairs,dtype=typ)
    
    dat_pairs[:,0] = dat_flt[pairs[:,0]]
    dat_pairs[:,1] = dat_flt[pairs[:,1]]
    
    return dat_pairs

def centroid_calc(mesh_size,origin,P,dz,hex_type):
    """
    
    Parameters
    ----------
    mesh size : iterable of int of size 3
        (nu,nv) - number of elements in the u and v directions
    origin : iterable of float of size 3
        (x0,y0,z0) - origin of the elements in the 
    P : flaot
        Flat-to-flat pitch of the hexagonal array
    dz : interable of float
        Gives the axial thickness of each layer in the hexagonal array
    hex_type : str
        Either 'x' or 'y' defining the orientation of the hexagonal array as 
        either x-type or y-type
    
    Returns
    -------
    uvw_arrs
        A set of arrays giving the u, v, and w location at each array element
    xyz_arrs
        A set of arrays giving the x, y, and, z centroid at each array element
    uvw_coords
        D
    xyz_coords
        D
    
    """
    # Get mesh size
    nu = mesh_size[0]
    nv = mesh_size[1]
    nw = mesh_size[2]
    
    # Hexagonal indices
    u = np.linspace(0,nu-1,nu,dtype=int).reshape(-1, 1, 1)
    v = np.linspace(0,nv-1,nv,dtype=int).reshape( 1,-1, 1)
    w = np.linspace(0,nw-1,nw,dtype=int).reshape( 1, 1,-1)
    
    # Centroids in cartesian coordinates as an array
    ua = 1*u + 0*v + 0*w
    va = 0*u + 1*v + 0*w
    wa = 0*u + 0*v + 1*w
    
    uvw_arrs = (ua,va,wa)
    
    # Hexagonal indices as a flat
    uf = ua.transpose([2,1,0]).reshape(-1,1)
    vf = va.transpose([2,1,0]).reshape(-1,1)
    wf = wa.transpose([2,1,0]).reshape(-1,1)
    uvw_flt = np.hstack([uf,vf,wf])
    
    # Projection of hexagonal indices onto cartesian coordinates
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
    
    # Centroids in cartesian coordinates as an array
    z = origin[2] + np.hstack([0, np.cumsum(dz[:-1]+dz[1:])/2 ])
    z = z.reshape(1,1,-1)
    
    xa = origin[0] + ux*u + vx*v + 0*z
    ya = origin[1] + uy*u + vy*v + 0*z
    za =              0*u +  0*v + 1*z
    
    xyz_arrs = (xa,ya,za)
    
    # Centroids as a flat
    xf = xa.transpose([2,1,0]).reshape(-1,1)
    yf = ya.transpose([2,1,0]).reshape(-1,1)
    zf = za.transpose([2,1,0]).reshape(-1,1)
    xyz_flt = np.hstack([xf,yf,zf])
    
    # Save centroids
    return uvw_arrs, xyz_arrs, uvw_flt, xyz_flt