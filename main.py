import numpy as np 
import sympy as sy
import tkinter as ttk
from tkinter.ttk import Combobox
from formulas import *

janela = ttk.Tk()
janela.title("Main")
janela.geometry("500x500")

def dimension():
  jnl = ttk.Tk()
  jnl.geometry("500x500")
  if cx == "1d":
    jnl.title("1D")

  if cx == "2d":
    jnl.title("2D")
    L_label = ttk.Label(jnl, text="L:", font=("Arial", 16))
    L_label.place(relx=0.1, rely=0.3, anchor=ttk.CENTER)
    L_entry = ttk.Entry(jnl, width=15, font=("Arial", 16))
    L_entry.place(relx=0.5, rely=0.3, anchor=ttk.E)

    nxi_label = ttk.Label(jnl, text="nxi:", font=("Arial", 16))
    nxi_label.place(relx=0.11, rely=0.4, anchor=ttk.CENTER)
    nxi_entry = ttk.Entry(jnl, width=2, font=("Arial", 16))
    nxi_entry.place(relx=0.21, rely=0.4, anchor=ttk.E)

    nyi_label = ttk.Label(jnl, text="nyi:", font=("Arial", 16))
    nyi_label.place(relx=0.11, rely=0.48, anchor=ttk.CENTER)
    nyi_entry = ttk.Entry(jnl, width=2, font=("Arial", 16))
    nyi_entry.place(relx=0.21, rely=0.48, anchor=ttk.E)

    nxf_label = ttk.Label(jnl, text="nxf:", font=("Arial", 16))
    nxf_label.place(relx=0.3, rely=0.4, anchor=ttk.CENTER)
    nxf_entry = ttk.Entry(jnl, width=2, font=("Arial", 16))
    nxf_entry.place(relx=0.4, rely=0.4, anchor=ttk.E)

    nyf_label = ttk.Label(jnl, text="nyf:", font=("Arial", 16))
    nyf_label.place(relx=0.3, rely=0.48, anchor=ttk.CENTER)
    nyf_entry = ttk.Entry(jnl, width=2, font=("Arial", 16))
    nyf_entry.place(relx=0.4, rely=0.48, anchor=ttk.E)

    def calc2d():
      L = float(L_entry.get())
      Nxi = int(nxi_entry.get())
      Nyi = int(nyi_entry.get())
      Nxf = int(nxf_entry.get())
      Nyf = int(nyf_entry.get())

      E = 2*Energia(hJ, m, L)/1.602E-19
      print(f"[A] - E: {np.format_float_scientific(E, precision = 2, exp_digits = 1)} eV\n")

      Eni = Energia_Inicial2D(Nxi, Nyi, hJ, m, L)

      print(f"[B] - Eni: {np.format_float_scientific(Eni/1.602E-19, precision = 2, exp_digits = 1)} eV\n")

      Enf = Energia_Final2D(Nxf, Nyf, hJ, m, L)

      print(f"[B] - Enf: {np.format_float_scientific(Enf/1.602E-19, precision = 2, exp_digits = 1)} eV\n")

      λi = Comprimento(hJ, m, Eni)
      λf = Comprimento(hJ, m, Enf)

      print(f"[C] - λi: {np.format_float_scientific(λi, precision = 2, exp_digits = 1)} m")
      print(f"[C] - λf: {np.format_float_scientific(λf, precision = 2, exp_digits = 1)} m\n")

    btn_calc = ttk.Button(jnl, text="Calcular", font=("Arial", 10), command=(calc2d))
    btn_calc.place(relx= 0.70, rely=0.4 , anchor=ttk.CENTER)


  if cx == "3d":
    jnl.title("3D")

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
  
def massa_eletron():
    global m
    m = Em
    selecionado_massa['text'] = "Massa: Eletron"

def massa_proton():
    global m
    m = Pm
    selecionado_massa['text'] = "Massa: Proton"

def def1d():
    global cx 
    cx = "1d"
    selecionado_d['text'] = "Dimesão: 1d"

def def2d():
    global cx 
    cx = "2d"
    selecionado_d['text'] = "Dimesão: 2d"

def def3d():
    global cx 
    cx = "3d"
    selecionado_d['text'] = "Dimesão: 3d"


titulo = ttk.Label(janela, text="Escolha a particula e a dimensão", font=("Arial Bold", 20))
titulo.place(relx=0.5, rely=0.2, anchor=ttk.CENTER)

massa_label = ttk.Label(janela, text="Massa:", font=("Arial Bold", 16))
massa_label.place(relx=0.15, rely=0.3, anchor=ttk.CENTER)

btn_eletron = ttk.Button(janela, text="Eletron", font=("Arial", 10), command=massa_eletron)
btn_eletron.place(relx= 0.30, rely=0.3 , anchor=ttk.CENTER)

btn_proton = ttk.Button(janela, text="Proton", font=("Arial", 10), command=massa_proton)
btn_proton.place(relx= 0.45, rely=0.3 , anchor=ttk.CENTER)

dimension_label = ttk.Label(janela, text="Dimensão da caixa:", font=("Arial Bold", 16))
dimension_label.place(relx=0.2655, rely=0.4, anchor=ttk.CENTER)

btn_1d = ttk.Button(janela, text="1D", font=("Arial", 10), command=def1d)
btn_1d.place(relx= 0.50, rely=0.4 , anchor=ttk.CENTER)

btn_2d = ttk.Button(janela, text="2D", font=("Arial", 10), command=def2d)
btn_2d.place(relx= 0.60, rely=0.4 , anchor=ttk.CENTER)

btn_3d = ttk.Button(janela, text="3D", font=("Arial", 10), command=def3d)
btn_3d.place(relx= 0.70, rely=0.4 , anchor=ttk.CENTER)

selecionado = ttk.Label(janela, text="Selecionados", font=("Arial", 10))
selecionado.place(relx=0.5, rely=0.75, anchor=ttk.CENTER)
selecionado_massa = ttk.Label(janela, text="Massa: ", font=("Arial", 10))
selecionado_massa.place(relx=0.5, rely=0.80, anchor=ttk.CENTER)
selecionado_d = ttk.Label(janela, text="Dimensão: ", font=("Arial", 10))
selecionado_d.place(relx=0.5, rely=0.85, anchor=ttk.CENTER)

continuar = ttk.Button(janela, text="Continuar", font=("Arial", 10), command=dimension)
continuar.place(relx= 0.5, rely=0.9, anchor=ttk.CENTER)

janela.mainloop()

print(m)
print(cx)

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
  
  
