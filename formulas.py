import numpy as np

def Energia(h, m, L):
    return (pow(h, 2))/(8*m*pow(L, 2))

def Energia_Inicial2D(nxi, nyi, h, m, L):
    return (pow(nxi, 2) + pow(nyi, 2))*Energia(h, m, L)

def Energia_Final2D(nxf, nyf, h, m, L):
    return (pow(nxf, 2) + pow(nyf, 2))*Energia(h, m, L)

def Energia_Inicial(ni, E):
    return (pow(ni, 2)*E)

def Comprimento(h, m, E):
    return h/np.sqrt(2*m*E)
