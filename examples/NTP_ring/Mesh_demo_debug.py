#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
from smeshtools.hextools import *
from smeshtools.hex_channel import Hex_channel
import matplotlib.pyplot as plt

# In[2]:


# Generate a text file that represents a hexagonal arrray with rings
channel_ids = ['F12','F11','F10','F09','F08','F07','R06','R05','R04','R03',
               'R02','R01']
ring_arr_2d = gen_hex_ring_2d_arr(channel_ids,void_id="VOI")
arr_to_hex_text(ring_arr_2d,"ring_hex_2d.txt",hex_type="x")


# In[3]:


# Verify that the text to arr reader returns the original array
ring_arr_2d_test = hex_text_to_arr("ring_hex_2d.txt",hex_type="x")
np.all(ring_arr_2d==ring_arr_2d_test)


# # 2. Create the channel-wise hex mesh model

# In[4]:


# Geometry specifications
P = 4.9
x0 = -P/4
y0 = -3*np.sqrt(3)/8*P
z0 = 0.5
origin = (x0,y0,z0)

# Axial geometry
z    = np.array([0, 5, 13, 21, 29, 37, 45, 53, 61, 69, 77, 85, 93, 101, 109, 117, 125, 135])
dz   = np.diff(z)
zmid = (z[1:]+z[:-1])/2


# In[16]:


# Input/Output Id map
layer_ids = np.linspace(1,17,17).astype(int).astype(str)
inp_id_map = np.zeros([0,17],dtype="U7")
for channel_id in channel_ids:
    channel_inp_ids = np.char.add(channel_id+'_',layer_ids)
    inp_id_map = np.vstack([inp_id_map, channel_inp_ids])


# In[24]:


hex_channel = Hex_channel("ring_hex_2d.txt", inp_id_map, channel_ids, 
                          out_id_map=None, hex_type="x", channel_void_id="VOI", 
                          inp_void_id="VOI", out_void_id="VOI", 
                          ele_void_id="VOI")
hex_channel.constrain(origin,P,dz)
hex_channel.gen_channel_centroids()

# Plot channel centroids
ch_centroids = hex_channel.channel_centroids
plt.figure()
for centroids in ch_centroids.values():
    plt.scatter(centroids[:,0],centroids[:,1])
plt.xlabel("x position")
plt.ylabel("y position")
plt.legend(ch_centroids.keys())
plt.show()

hex_channel.gen_channel_GeometryCollection()

# In[]:

hex_channel.gen_interfaces()