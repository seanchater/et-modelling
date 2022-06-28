import numpy as np

vc = ""
ndvi = ""


# **effective_leaf_area_index**************************************************
# constants or predefined:
""" nd_min = np.nanmin(ndvi)
nd_max = np.nanmax(ndvi)
vc_pow = 0.7
vc_min = np.nanmin(vc)
vc_max = np.nanmax(vc)
lai_pow = -0.45 """



def exception(test):
    try:
        float(test)
        return True
    except:
        return False

lai_pow = 0.0
vc_pow = 0.0

def setPowLAI(value):
    global lai_pow
    lai_pow = value

def getPowLAI():
    global lai_pow
    return lai_pow

in_ = input("LAI pow: ")
setPowLAI([lambda:-0.45, lambda: float(in_)][exception(in_)]())
print(getPowLAI())
print(exception(in_))