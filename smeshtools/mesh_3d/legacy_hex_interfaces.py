# -*- coding: utf-8 -*-
"""
Created on Thu May 15 00:51:14 2025

@author: 17066
"""

from collections import Counter
import numpy as np

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