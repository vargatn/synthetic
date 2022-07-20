

import numpy as np



def scale_image(canvas):
    try:
        res = np.arcsinh(canvas) / canvas
    except:
        res = np.arcsinh(canvas.array) / canvas.array
    return