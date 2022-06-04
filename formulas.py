import numpy as np

def Energia(h, m, L):
    return (pow(h, 2))/(8*m*pow(L, 2))

def Energia_Inicial(ni, E):
    return (pow(ni, 2)*E)

def Energia_Final(nf, E):
    return (pow(nf, 2)*E)

def Energia_Inicial2D(nxi, nyi, h, m, L):
    return (pow(nxi, 2) + pow(nyi, 2))*Energia(h, m, L)

def Energia_Final2D(nxf, nyf, h, m, L):
    return (pow(nxf, 2) + pow(nyf, 2))*Energia(h, m, L)

def Energia_Inicial3D(nxi, nyi, nzi, h, m, L):
  return (pow(nxi, 2) + pow(nyi, 2) + pow(nzi, 2))*Energia(h, m, L)

def Energia_Final3D(nxf, nyf, nzf, h, m, L):
  return (pow(nxf, 2) + pow(nyf, 2) + pow(nzf, 2))*Energia(h, m, L)

def Comprimento(h, m, E):
    return h/np.sqrt(2*m*E)

def Comprimento_Foton(h, E):
  return (h*3e8)/E

def A_calc(L):
  return np.sqrt(2/L)

def L_calc(A):
  return 2/pow(A, 2)

def n_calc(k, L):
  return (k*L)/np.pi

def k_Inicial(ni, L):
  return (ni*np.pi)/L

def k_Final(nf, L):
  return (nf*np.pi)/L

def Velocidade(E, m):
  return np.sqrt((2*E)/m)

def Theta(ni, a, L):
  return ni*np.pi*a/L
