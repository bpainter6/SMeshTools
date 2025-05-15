# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 16:54:39 2024

@author: 17066
"""

import numpy as np
from copy import deepcopy
from textwrap import wrap

def gen_hex_2d_arr(nu,nv):
    """Generate an array representing data on a hexagonal grid. Each element
    is a string representing a unique number.
    
    Parameters
    ----------
    nu - int
        number of elements along axis 0
    nv - int
        number of elements along axis 1
    
    """
    
    # Number of elements
    n = nu*nv
    
    # Generate the array
    arr = np.linspace(1,n,n).astype(int)
    arr = np.vectorize(str)(arr)
    arr = arr.reshape(nv,nu).transpose([1,0])
    
    return arr

def gen_hex_3d_arr(nu,nv,nw):
    """Generate an array representing data on a hexagonal grid. Each element
    is a string representing a unique number.
    
    Parameters
    ----------
    nu - int
        number of elements along axis 0
    nv - int
        number of elements along axis 1
    nw - int
        number of elements along axis 2
    
    """
    
    # Number of elements
    n = nu*nv*nw
    
    # Generate the array
    arr = np.linspace(1,n,n).astype(int)
    arr = np.vectorize(str)(arr)
    arr = arr.reshape(nw,nv,nu).transpose([2,1,0])
    
    return arr

def gen_hex_ring_2d_arr(ring_ids,void_id="VD"):
    """Generate an array representing a hexagonal geometry where each ring 
    is given the same id provided by ring_ids (from center to out)
    
    Parameters
    ----------
    ring_ids - 1d sequence of str
        ids given to each ring from center to out
    void_id - str
        id given to regions associated with a void
    
    """
    
    nr = len(ring_ids)
    
    # Calculate the upper left and bottom right components of the array
    A = np.full([nr,nr],"",dtype='<U7')
    for i,ring_id in enumerate(ring_ids):
        vec = np.array([ring_id]*(i+1))
        A = np.char.add(A,np.diag(vec,nr-i-1))
    A = A[::-1,:]
    UL = deepcopy(A)
    BR = A[::-1,::-1][1:,1:]
    
    # Calculate the upper right and bottom left components
    B = np.full([nr,nr],"",dtype='<U7')
    for i,ring_id in enumerate(ring_ids[::-1]):
        B[:nr-i,:nr-i] = ring_id
    B = B[::-1,:]
    UR = B[:,1:]
    BL = B[::-1,::-1][1:,:]
    
    # Calculate the block array
    arr = np.block([[UL,UR],
                    [BL,BR]])
    
    # Insert void ids
    mask = np.where(arr=="")
    arr[mask] = void_id
    
    return arr

def gen_hex_ring_3d_arr(ring_ids,void_id="VD"):
    """Generate an array representing a hexagonal geometry where each ring 
    is given the same id provided by ring_ids (from center to out)
    
    Parameters
    ----------
    ring_ids - 2d sequence of str
        ids given to each ring from center to out
    void_id - str
        id given to regions associated with a void
    
    """
    
    ring_ids = np.array(ring_ids)
    nz,nr = np.shape(ring_ids)
    arr = np.full([2*nr-1,2*nr-1,nz]," ",dtype='<U7')
    
    # create ring ids for each layer
    for i in range(nz):
        arr[:,:,i] = gen_hex_ring_2d_arr(ring_ids[i,:],void_id)
    
    return arr

def arr_pad_to_arr_txt(arr_pad,space_len):
    ny,nx = np.shape(arr_pad)
    arr_txt = np.full([ny,2*nx+ny-2]," "*space_len,dtype='<U7')
    for i in range(ny):
        lb = i
        ub = 2*nx+i
        arr_txt[i,lb:ub:2] = arr_pad[i,:]
    return arr_txt
    

def arr_to_hex_text(arr,filename,hex_type="x",txt_len_min="auto",pad_val=" ",
                    space_len="auto"):
    """Converts an array representation of a hexagonal mesh to a text file 
    where it is displayed as a hexagonal lattice
    
    Parameters
    ---------
    arr - numpy array
        array representing the 
    filename - str
        path to write the hex text
    hex_type - str
        specifies oriendtation of data printout
    txt_len_min - int/str/None
        minimum length of the written data. Padding is applied to makeup any
        difference
        "auto" - automatically calculates the minimum padding necessary
        1 - disables padding
    pad_val - txt
        string used to pad written data
    space_len - int/str
        length of the space between written data
        "auto" - automatically makes the spacing equal to the minimum text
        length
    
    """
    
    # Legnth of data when converted to str
    def data_len(item):
        return len(str(item))
    
    # Calculate minimum text length
    if txt_len_min=="auto":
        txt_len_min = np.max(np.vectorize(data_len)(arr))
    
    # Calculate space length
    if space_len=="auto":
        space_len = deepcopy(txt_len_min)
    
    # Padding function
    def pad(item):
        item = str(item)
        left = max(txt_len_min-len(item),0)
        return pad_val*left+item
    
    # Pad text
    arr_pad = np.vectorize(pad)(arr)
    
    # Calculate shape
    ndim = arr_pad.ndim
    if ndim==2:
        nu,nv = np.shape(arr_pad)
    
    elif ndim==3:
        nu,nv,nw = np.shape(arr_pad)
    
    # Reorient array depending on type
    if ndim==2:
        if hex_type=="x":
            arr_pad = arr_pad.transpose([1,0])
    
    elif ndim==3:
        if hex_type=="x":
            arr_pad = arr_pad.transpose([1,0,2])
    
    # Format user provided array as a checkered numpy array
    arr_txts = []
    if ndim==2:
        arr_txts.append(arr_pad_to_arr_txt(arr_pad,space_len))
    elif ndim==3:
        for i in range(arr_pad.shape[-1]):
            arr_padi = arr_pad[:,:,i]
            arr_txts.append(arr_pad_to_arr_txt(arr_padi,space_len))
    
    # Print results to file
    f = open(filename, "w")
    for arr_txt in arr_txts:
        for row in arr_txt[::-1]:
            line = ''.join(row)
            line += "\n"
            f.write(line)
        f.write("\n")
    f.close()

def arr_text_to_arr(arr_text,hex_type):
    arr = np.array(arr_text)
    ny,nx = np.shape(arr)
    arr = arr[::-1]
    
    # reshape data
    if hex_type=="x":
        arr = arr.transpose([1,0])
    
    return arr

def hex_text_to_arr(filename,hex_type="x",ndim="auto"):
    """Converts an array representation of a hexagonal mesh to a text file 
    where it is displayed in a x-hexagonal lattice"""
    
    arr_text = []
    arrs = []
    f = open(filename, "r")
    
    for line in f.readlines():
        if line=="\n":
            arr = arr_text_to_arr(arr_text,hex_type)
            arrs.append(arr)
            arr_text = []
        else:
            arr_text.append(line.split())
    
    # Post process that last array
    if len(arr_text)>0:
        arr = arr_text_to_arr(arr_text,hex_type)
        arrs.append(arr)
    f.close()
    
    if (ndim=="auto" and len(arrs)==1) or ndim==2:
        return arrs[0]
    elif (ndim=="auto" and len(arrs)>1) or ndim==3:
        arr = np.array(arrs)
        return arr.transpose([1,2,0])