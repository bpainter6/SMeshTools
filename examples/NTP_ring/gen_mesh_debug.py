# -*- coding: utf-8 -*-
"""
Created on Mon Oct 28 17:56:41 2024

@author: 17066
"""

from xsfit.preprocessing.hextools import gen_arr,arr_to_y_hex_text,hex_text_to_arr
arr = gen_arr(5,True,2)
arr_to_y_hex_text(arr,"hex_mesh.txt",2)
arr2 = hex_text_to_arr("hex_mesh.txt",2)