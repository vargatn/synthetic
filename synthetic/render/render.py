"""
Galsim rendering backend

the first stretch goal is to build a stand alone renderer which is encapsulated in a class. The 

# TODO add input catalog specification here

# TODO add color composite image creator function



"""

import numpy as np



def scale_image(canvas):
    try:
        res = np.arcsinh(canvas) / canvas
    except:
        res = np.arcsinh(canvas.array) / canvas.array
    return






def color_composite():
    pass

