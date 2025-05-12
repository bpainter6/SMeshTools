# -*- coding: utf-8 -*-
"""
Created on Wed Apr  2 23:29:49 2025

@author: 17066
"""

import numpy as np
from smeshtools.hextools import *

ring_ids_2d = ["a","b","c","d","e","f","g"]
ring_arr_2d = gen_hex_ring_2d_arr(ring_ids_2d)
arr_to_hex_text(ring_arr_2d,"ring_hex_2d.txt",hex_type="x")

x_arr_2d = gen_hex_2d_arr(3,6)
arr_to_hex_text(x_arr_2d,"x_hex_2d.txt",hex_type="x")
x_arr_2d_test = hex_text_to_arr("x_hex_2d.txt",hex_type="x")

y_arr_2d = gen_hex_2d_arr(8,2)
arr_to_hex_text(y_arr_2d,"y_hex_2d.txt",hex_type="y")
y_arr_2d_test = hex_text_to_arr("y_hex_2d.txt",hex_type="y")

ring_ids_3d = [["a_1","b_1","c_1","d_1","e_1","f_1","g_1"],
               ["a_2","b_2","c_2","d_2","e_2","f_2","g_2"],
               ["a_3","b_3","c_3","d_3","e_3","f_3","g_3"]]
ring_arr_3d = gen_hex_ring_3d_arr(ring_ids_3d)
arr_to_hex_text(ring_arr_3d,"ring_hex_3d.txt",hex_type="x")
ring_arr_3d_test = hex_text_to_arr("ring_hex_3d.txt",hex_type="x")

x_arr_3d = gen_hex_3d_arr(8,2,5)
arr_to_hex_text(x_arr_3d,"x_hex_3d.txt",hex_type="x")
x_arr_3d_test = hex_text_to_arr("x_hex_3d.txt",hex_type="x")

y_arr_3d = gen_hex_3d_arr(3,6,3)
arr_to_hex_text(y_arr_3d,"y_hex_3d.txt",hex_type="y")
y_arr_3d_test = hex_text_to_arr("y_hex_3d.txt",hex_type="y")