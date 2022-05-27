import numpy as np 
import sympy as sy

def fi(x):
  return (sy.abs(sy.sqrt(2/L) * sy.sin((ni*sy.pi()/L)*x)))**2
 
def ff(x):
  return (sy.abs(sy.sqrt(2/L) * sy.sin((nf*np.pi()/L)*x)))**2

m = 0 # Massa
Em = 9.10e-31 # Massa Eletron
Pm = 1.67e-27 # Massa Proton
hJ = 6.626e-34 # Constante de Planck em Joule
hEv = 4.136e-15 # Constante de Planck em Ev
c = 3e8 # Velocidade da Luz

type = input("Particula confinada: [Eletron|Proton]")
type = type.lower()

if type == "eletron":
  m = Em
if type == "proton":
  m = Pm

cx = input("Tipo de caixa: [1D|2D|3D]")
cx = cx.lower()

if cx == "1d":
  print("[1] - Determinação da função de onda quântica e outros parâmetros")
  print("[2] - Cálculo dos parâmetros da caixa e partícula, dada a função de onda")
  t = int(input())
  if t == 1:
    print("Se não houver valor, digite null")
    L = float(input("Largura: "))
    ni = int(input("Digite o número quantico inicial: "))
    nf = int(input("Digite o número quantico final: "))

    E = (pow(hJ, 2)/((8*m) * pow(L, 2)))
    print(f"[A] E1: {np.format_float_scientific(E, precision = 2, exp_digits = 1)}")
    print("Região do poço para procurar pela partícula [para cálculo de P(a ≤ x ≤ b)]: ")

    a = float(input("a - Coordenada dentro do poço: "))

    b = float(input("b - Coordenada dentro do poço: "))  

    x = sy.Symbol("x")

    A = np.sqrt(2/L)

    Ki = (ni*np.pi)/L
    Kf = (nf*np.pi)/L

    print(f"A: {np.format_float_scientific(A, precision=3, exp_digits=1)} Ki: {np.format_float_scientific(Ki, precision=3, exp_digits=1)} Kf {np.format_float_scientific(Kf, precision=3, exp_digits=1)}\n")

    # ----- # A # ----- #

    def Func_quant_ini(x):
     return A*np.sin(Ki*np.pi*ni*x)
    
    def Func_quant_fin(x):
     return A*np.sin(Kf*np.pi*nf*x)

    # ----- # B # ----- #

    Eni = (pow(ni, 2)*E)
    Enf = (pow(nf, 2)*E)
    
    print(f"[B] - Eni: {np.format_float_scientific(Eni, precision=3, exp_digits=1)} J | Enf: {np.format_float_scientific(Enf, precision=3, exp_digits=1)} J")

    print(f"[B] - Eni: {np.format_float_scientific(Eni/1.602e-19, precision=3, exp_digits=1)} Ev | Enf: {np.format_float_scientific(Enf/1.602e-19, precision=3, exp_digits=1)} Ev\n")
    
    
    En = Enf/1.602e-19 - Eni/1.602e-19
    λ = hEv*c/En

    print(f"[C] - En: {(np.format_float_scientific(En, precision=3, exp_digits=1))} Ev | En: {(np.format_float_scientific(λ, precision=3, exp_digits=1))} J\n") #C

    Vi = np.sqrt((2*Eni)/m)
    Vf = np.sqrt((2*Enf)/m)

    #D
    print(f"[D] - Vi: {(np.format_float_scientific(Vi, precision=3, exp_digits=1))} | Vf: {(np.format_float_scientific(Vf, precision=3, exp_digits=1))}\n")

      # ----- # E # ----- #

    λi = hJ/np.sqrt(2*m*Eni)
    λf = hJ/np.sqrt(2*m*Enf)

    print(f"[E] - λi: {(np.format_float_scientific(λi, precision=3, exp_digits=1))} | λf: {(np.format_float_scientific(λf, precision=3, exp_digits=1))}\n")

     # ----- # F # ----- #


    thetai = ni*np.pi*a/L
    thetaf = ni*np.pi*b/L
    p = sy.integrate((2/ni*sy.pi)*(pow(sy.sin(x),2)), (x, thetai, thetaf))
    print(f"[F] - ni = {np.format_float_scientific(ni, precision = 2, exp_digits = 1)}, P({np.format_float_scientific(thetai, precision = 2, exp_digits = 1)} ≤ x ≤ {np.format_float_scientific(thetaf, precision = 2, exp_digits = 1)}) = {np.format_float_scientific(p*10, precision = 2, exp_digits = 1)}%\n")
    thetai = nf*np.pi*a/L
    thetaf = nf*np.pi*b/L
    p = sy.integrate((2/nf*sy.pi)*(pow(sy.sin(x),2)), (x, thetai, thetaf))
    print(f"[F] - nf = {np.format_float_scientific(nf, precision = 2, exp_digits = 1)}, P({np.format_float_scientific(thetai, precision = 2, exp_digits = 1)} ≤ x ≤ {np.format_float_scientific(thetaf, precision = 2, exp_digits = 1)}) = {np.format_float_scientific(p*10, precision = 2, exp_digits = 1)}%")
  if t == 2:
    A = float(input("Digite um valor de A no SI: "))
    k = float(input("Digite um valor de k no SI: "))

    # 2 ----- # A # ----- #

    L = 2/pow(A, 2)

    print(f"[A] - Largura: {np.format_float_scientific(L, precision = 2, exp_digits = 1)} m | {np.format_float_scientific(L*10e8, precision = 2, exp_digits = 1)} nm\n")

    # 2 ----- # B # ----- #

    n = (k*L)/np.pi

    print(f"[B] - n: {round(n)}\n")

    # 2 ----- # C # ----- #

    E = (pow(hJ, 2)/((8*m) * pow(L, 2)))
    

    print(f"[C] - E: {np.format_float_scientific(E, precision = 2, exp_digits = 1)} J | E: {np.format_float_scientific(E/1.602e-19, precision = 2, exp_digits = 1)} Ev")

    # 2 ----- # D # ----- #

    V = np.sqrt((2*E)/m)

    print(f"[E] - V: {np.format_float_scientific(V, precision = 2, exp_digits = 1)}\n")

    #  ----- # Exercicio 3 # ----- #
if cx == "2d":
  L = float(input("Largura: "))

  Nxi = int(input("Numero inicial da particula confinada x: "))
  Nyi = int(input("Numero inicial da particula confinada y: "))

  Nxf = int(input("Numero final da particula confinada x: "))
  Nyf = int(input("Numero final da particula confinada y: "))

  # 3 ----- # A # ----- #

  E = ((2)*(pow(hJ, 2))/(8*m*pow(L, 2)))/1.602E-19
  print(f"[A] - E: {np.format_float_scientific(E, precision = 2, exp_digits = 1)} eV\n")

  Eni = ((pow(Nxi, 2) + pow(Nyi, 2))*(pow(hJ, 2))/(8*m*pow(L, 2)))

  print(f"[B] - Eni: {np.format_float_scientific(Eni/1.602E-19, precision = 2, exp_digits = 1)} eV\n")

  Enf = ((pow(Nxf, 2) + pow(Nyf, 2))*(pow(hJ, 2))/(8*m*pow(L, 2)))

  print(f"[B] - Enf: {np.format_float_scientific(Enf/1.602E-19, precision = 2, exp_digits = 1)} eV\n")

  λi = hJ/np.sqrt(2*m*Eni)
  λf = hJ/np.sqrt(2*m*Enf)

  print(f"[C] - λi: {np.format_float_scientific(λi, precision = 2, exp_digits = 1)} m")
  print(f"[C] - λf: {np.format_float_scientific(λf, precision = 2, exp_digits = 1)} m\n")

  


 #  ----- # Exercicio 4 # ----- #
if cx == "3d":
  L = float(input("Largura: "))

  Nxi = int(input("Numero inicial da particula confinada x: "))
  Nyi = int(input("Numero inicial da particula confinada y: "))
  Nzi = int(input("Numero inicial da particula confinada z: "))

  Nxf = int(input("Numero final da particula confinada x: "))
  Nyf = int(input("Numero final da particula confinada y: "))
  Nzf = int(input("Numero final da particula confinada z: "))

  E = ((3)*(pow(hJ,2))/(8*m*pow(L, 2)))/1.602E-19
  print(f"[A] - {np.format_float_scientific(E, precision = 2, exp_digits = 1)}")

  Eni = ((pow(Nxi, 2) + pow(Nyi, 2) + pow(Nzi, 2))*(pow(hJ, 2))/(8*m*pow(L, 2)))/1.602E-19

  print(f"[B] - Eni: {np.format_float_scientific(Eni, precision = 2, exp_digits = 1)} eV")

  Enf = ((pow(Nxf, 2) + pow(Nyf, 2) + pow(Nzf, 2))*(pow(hJ, 2))/(8*m*pow(L, 2)))/1.602E-19

  print(f"[B] - Enf: {np.format_float_scientific(Enf, precision = 2, exp_digits = 1)} eV")

  En = Eni - Enf

  λ = (hEv*c)/En

  print(f"[C] - En: {np.format_float_scientific(En, precision = 2, exp_digits = 1)} eV")
  
  print(f"[C] - λ: {np.format_float_scientific(λ, precision = 2, exp_digits = 1)} eV")
  
  
