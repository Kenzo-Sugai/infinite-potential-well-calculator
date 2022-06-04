import numpy as np 
import sympy as sy
import tkinter as ttk
from tkinter.ttk import Combobox
from formulas import *

m = 0 # Massa
Em = 9.10e-31 # Massa Eletron
Pm = 1.67e-27 # Massa Proton
hJ = 6.626e-34 # Constante de Planck em Joule
hEv = 4.136e-15 # Constante de Planck em Ev
c = 3e8 # Velocidade da Luz

janela = ttk.Tk()
janela.title("Main")
janela.geometry("500x500")

def dimension():
  jnl = ttk.Tk()
  jnl.geometry("500x500")
  if cx == "1d":
    jnl.title("1D")
    L_label = ttk.Label(jnl, text="L:", font=("Arial", 16))
    L_label.place(relx=0.1, rely=0.2, anchor=ttk.CENTER)
    L_entry = ttk.Entry(jnl, width=15, font=("Arial", 16))
    L_entry.place(relx=0.5, rely=0.2, anchor=ttk.E)

    ni_label = ttk.Label(jnl, text="ni:", font=("Arial", 16))
    ni_label.place(relx=0.11, rely=0.3, anchor=ttk.CENTER)
    ni_entry = ttk.Entry(jnl, width=2, font=("Arial", 16))
    ni_entry.place(relx=0.21, rely=0.3, anchor=ttk.E)

    nf_label = ttk.Label(jnl, text="nf:", font=("Arial", 16))
    nf_label.place(relx=0.11, rely=0.38, anchor=ttk.CENTER)
    nf_entry = ttk.Entry(jnl, width=2, font=("Arial", 16))
    nf_entry.place(relx=0.21, rely=0.38, anchor=ttk.E)

    a_label = ttk.Label(jnl, text="a:", font=("Arial", 16))
    a_label.place(relx=0.25, rely=0.3, anchor=ttk.CENTER)
    a_entry = ttk.Entry(jnl, width=8, font=("Arial", 16))
    a_entry.place(relx=0.50, rely=0.3, anchor=ttk.E)

    b_label = ttk.Label(jnl, text="b:", font=("Arial", 16))
    b_label.place(relx=0.25, rely=0.38, anchor=ttk.CENTER)
    b_entry = ttk.Entry(jnl, width=8, font=("Arial", 16))
    b_entry.place(relx=0.50, rely=0.38, anchor=ttk.E)

    A_label = ttk.Label(jnl, text="A:", font=("Arial", 16))
    A_label.place(relx=0.1, rely=0.6, anchor=ttk.CENTER)
    A_entry = ttk.Entry(jnl, width=15, font=("Arial", 16))
    A_entry.place(relx=0.5, rely=0.6, anchor=ttk.E)

    k_label = ttk.Label(jnl, text="k:", font=("Arial", 16))
    k_label.place(relx=0.1, rely=0.68, anchor=ttk.CENTER)
    k_entry = ttk.Entry(jnl, width=15, font=("Arial", 16))
    k_entry.place(relx=0.5, rely=0.68, anchor=ttk.E)

    def calc1dA():
      L = float(L_entry.get())
      ni = int(ni_entry.get())
      nf = int(nf_entry.get())

      E = Energia(hJ, m, L)

      a = float(a_entry.get())

      b = float(b_entry.get())  

      x = sy.Symbol("x")

      A = A_calc(L)

      Ki = k_Inicial(ni, L)
      Kf = k_Final(nf, L)

      K = Kf - Ki

      # ----- # A # ----- #

      def Func_quant_ini(x):
        return A*np.sin(Ki*np.pi*ni*x)
      
      def Func_quant_fin(x):
        return A*np.sin(Kf*np.pi*nf*x)

      # ----- # B # ----- #

      Eni = Energia_Inicial(ni, E)
      Enf = Energia_Final(nf, E)
      
      En = Enf/1.602e-19 - Eni/1.602e-19

      λ = Comprimento_Foton(hEv, En)

      Vi = Velocidade(Eni, m)
      Vf = Velocidade(Enf, m)

      #D

        # ----- # E # ----- #

      λi = Comprimento(hJ, m, Eni)
      λf = Comprimento(hJ, m, Enf)

      # ----- # F # ----- #

      thetai = Theta(ni, a, L)
      thetaf = Theta(ni, b, L)

      pi = sy.integrate((2/ni*sy.pi)*(pow(sy.sin(x),2)), (x, thetai, thetaf))
      
      thetai = Theta(nf, a, L)
      thetaf = Theta(nf, b, L)
      pf = sy.integrate((2/nf*sy.pi)*(pow(sy.sin(x),2)), (x, thetai, thetaf))

      resp = ttk.Tk()
      resp.geometry("500x500")
      resp.title("Resposta")
      
      Energia_label = ttk.Label(resp, text=f"Energia Fundamental: {np.format_float_scientific(E, precision = 2, exp_digits = 1)}", font=("Arial", 16))
      Energia_label.place(relx=0.5, rely=0.15, anchor=ttk.CENTER)

      form_label = ttk.Label(resp, text=f"A: {np.format_float_scientific(A, precision=3, exp_digits=1)} Ki: {np.format_float_scientific(Ki, precision=3, exp_digits=1)} Kf: {np.format_float_scientific(Kf, precision=3, exp_digits=1)}", font=("Arial", 16))
      form_label.place(relx=0.5, rely=0.20, anchor=ttk.CENTER)

      Energia_Inicial_label = ttk.Label(resp, text=f"Energia Incial: {np.format_float_scientific(Eni, precision=3, exp_digits=1)} J | {np.format_float_scientific(Eni/1.602E-19, precision = 2, exp_digits = 1)} Ev", font=("Arial", 16))
      Energia_Inicial_label.place(relx=0.5, rely=0.25, anchor=ttk.CENTER)

      Energia_Final_label = ttk.Label(resp, text=f"Energia Final: {np.format_float_scientific(Enf, precision=3, exp_digits=1)} J | {np.format_float_scientific(Enf/1.602E-19, precision = 2, exp_digits = 1)} Ev", font=("Arial", 16))
      Energia_Final_label.place(relx=0.5, rely=0.30, anchor=ttk.CENTER)

      Energia_Foton_label = ttk.Label(resp, text=f"Energia do Fóton: {np.format_float_scientific(En, precision = 2, exp_digits = 1)} Ev", font=("Arial", 16))
      Energia_Foton_label.place(relx=0.5, rely=0.35, anchor=ttk.CENTER)

      Comprimento_Foton_label = ttk.Label(resp, text=f"Comprimento do Fóton: {np.format_float_scientific(λ, precision = 2, exp_digits = 1)} m", font=("Arial", 16))
      Comprimento_Foton_label.place(relx=0.5, rely=0.40, anchor=ttk.CENTER)

      Velocidade_Inicial_label = ttk.Label(resp, text=f"Velocidade Inicial: {np.format_float_scientific(Vi, precision = 2, exp_digits = 1)} m/s", font=("Arial", 16))
      Velocidade_Inicial_label.place(relx=0.5, rely=0.45, anchor=ttk.CENTER)

      Velocidade_Final_label = ttk.Label(resp, text=f"Velocidade Final: {np.format_float_scientific(Vf, precision = 2, exp_digits = 1)} m/s", font=("Arial", 16))
      Velocidade_Final_label.place(relx=0.5, rely=0.50, anchor=ttk.CENTER)

      Comprimento_Inicial_label = ttk.Label(resp, text=f"Comprimento Inicial: {np.format_float_scientific(λi, precision = 2, exp_digits = 1)} m", font=("Arial", 16))
      Comprimento_Inicial_label.place(relx=0.5, rely=0.55, anchor=ttk.CENTER)

      Comprimento_Final_label = ttk.Label(resp, text=f"Comprimento Final: {np.format_float_scientific(λf, precision = 2, exp_digits = 1)} m", font=("Arial", 16))
      Comprimento_Final_label.place(relx=0.5, rely=0.60, anchor=ttk.CENTER)

      Probabilidade_Inicial_label = ttk.Label(resp, text=f"Probabilidade Inicial: {np.format_float_scientific(pi, precision = 2, exp_digits = 1)} %", font=("Arial", 16))
      Probabilidade_Inicial_label.place(relx=0.5, rely=0.65, anchor=ttk.CENTER)

      Probabilidade_Final_label = ttk.Label(resp, text=f"Probabilidade Final: {np.format_float_scientific(pf, precision = 2, exp_digits = 1)} %", font=("Arial", 16))
      Probabilidade_Final_label.place(relx=0.5, rely=0.70, anchor=ttk.CENTER)

    def calc1dB():
      A = float(A_entry.get())
      k = float(k_entry.get())

      # 2 ----- # A # ----- #

      L = L_calc(A)

      # 2 ----- # B # ----- #

      n = n_calc(k, L)

      # 2 ----- # C # ----- #

      E = Energia(hJ, m, L)

      # 2 ----- # D # ----- #

      V = Velocidade(E, m)

      resp = ttk.Tk()
      resp.geometry("500x500")
      resp.title("Resposta")

      Largura_label = ttk.Label(resp, text=f"Largura: {np.format_float_scientific(L*10e8, precision=3, exp_digits=1)} nm", font=("Arial", 16))
      Largura_label.place(relx=0.5, rely=0.10, anchor=ttk.CENTER)

      n_label = ttk.Label(resp, text=f"n: {np.format_float_scientific(n, precision=3, exp_digits=1)}", font=("Arial", 16))
      n_label.place(relx=0.5, rely=0.20, anchor=ttk.CENTER)

      Energia_label = ttk.Label(resp, text=f"Energia Fundamental: {np.format_float_scientific(E, precision = 2, exp_digits = 1)} J | {np.format_float_scientific(E/1.602e-19, precision = 2, exp_digits = 1)} Ev", font=("Arial", 16))
      Energia_label.place(relx=0.5, rely=0.30, anchor=ttk.CENTER)

      Velocidade_label = ttk.Label(resp, text=f"Energia Final: {np.format_float_scientific(V, precision=3, exp_digits=1)} m/s", font=("Arial", 16))
      Velocidade_label.place(relx=0.5, rely=0.40, anchor=ttk.CENTER)

    btn_calcA = ttk.Button(jnl, text="Calcular", font=("Arial", 10), command=(calc1dA))
    btn_calcA.place(relx= 0.65, rely=0.3 , anchor=ttk.CENTER)

    btn_calcB = ttk.Button(jnl, text="Calcular", font=("Arial", 10), command=(calc1dB))
    btn_calcB.place(relx= 0.65, rely=0.6 , anchor=ttk.CENTER)
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

      Eni = Energia_Inicial2D(Nxi, Nyi, hJ, m, L)

      Enf = Energia_Final2D(Nxf, Nyf, hJ, m, L)

      λi = Comprimento(hJ, m, Eni)
      λf = Comprimento(hJ, m, Enf)

      resp = ttk.Tk()
      resp.geometry("500x500")
      resp.title("Resposta")

      Energia_label = ttk.Label(resp, text=f"Energia Fundamental: {np.format_float_scientific(E, precision = 2, exp_digits = 1)}", font=("Arial", 16))
      Energia_label.place(relx=0.5, rely=0.1, anchor=ttk.CENTER)

      Energia_Inicial_label = ttk.Label(resp, text=f"Energia Incial: {np.format_float_scientific(Eni/1.602E-19, precision = 2, exp_digits = 1)} Ev", font=("Arial", 16))
      Energia_Inicial_label.place(relx=0.5, rely=0.3, anchor=ttk.CENTER)

      Energia_Final_label = ttk.Label(resp, text=f"Energia Final: {np.format_float_scientific(Enf/1.602E-19, precision = 2, exp_digits = 1)} Ev", font=("Arial", 16))
      Energia_Final_label.place(relx=0.5, rely=0.5, anchor=ttk.CENTER)

      Comprimento_Inicial_label = ttk.Label(resp, text=f"Comprimento Incial: {np.format_float_scientific(λi, precision = 2, exp_digits = 1)} m", font=("Arial", 16))
      Comprimento_Inicial_label.place(relx=0.5, rely=0.7, anchor=ttk.CENTER)

      Comprimento_Final_label = ttk.Label(resp, text=f"Comprimento Final: {np.format_float_scientific(λf, precision = 2, exp_digits = 1)} m", font=("Arial", 16))
      Comprimento_Final_label.place(relx=0.5, rely=0.9, anchor=ttk.CENTER)

    btn_calc = ttk.Button(jnl, text="Calcular", font=("Arial", 10), command=(calc2d))
    btn_calc.place(relx= 0.70, rely=0.4 , anchor=ttk.CENTER)


  if cx == "3d":
    jnl.title("3D")
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

    nzi_label = ttk.Label(jnl, text="nzi:", font=("Arial", 16))
    nzi_label.place(relx=0.11, rely=0.56, anchor=ttk.CENTER)
    nzi_entry = ttk.Entry(jnl, width=2, font=("Arial", 16))
    nzi_entry.place(relx=0.21, rely=0.56, anchor=ttk.E)

    nxf_label = ttk.Label(jnl, text="nxf:", font=("Arial", 16))
    nxf_label.place(relx=0.3, rely=0.4, anchor=ttk.CENTER)
    nxf_entry = ttk.Entry(jnl, width=2, font=("Arial", 16))
    nxf_entry.place(relx=0.4, rely=0.4, anchor=ttk.E)

    nyf_label = ttk.Label(jnl, text="nyf:", font=("Arial", 16))
    nyf_label.place(relx=0.3, rely=0.48, anchor=ttk.CENTER)
    nyf_entry = ttk.Entry(jnl, width=2, font=("Arial", 16))
    nyf_entry.place(relx=0.4, rely=0.48, anchor=ttk.E)

    nzf_label = ttk.Label(jnl, text="nzf:", font=("Arial", 16))
    nzf_label.place(relx=0.3, rely=0.56, anchor=ttk.CENTER)
    nzf_entry = ttk.Entry(jnl, width=2, font=("Arial", 16))
    nzf_entry.place(relx=0.4, rely=0.56, anchor=ttk.E)

    def calc3d():
      L = float(L_entry.get())
      Nxi = int(nxi_entry.get())
      Nyi = int(nyi_entry.get())
      Nzi = int(nzi_entry.get())
      Nxf = int(nxf_entry.get())
      Nyf = int(nyf_entry.get())
      Nzf = int(nzf_entry.get())

      E = 3*Energia(hJ, m, L)/1.602E-19

      Eni = Energia_Inicial3D(Nxi, Nyi, Nzi, hJ, m, L)/1.602E-19

      Enf = Energia_Final3D(Nxf, Nyf, Nzf, hJ, m, L)/1.602E-19

      En = Eni - Enf

      λ = Comprimento_Foton(hEv, En)

      resp = ttk.Tk()
      resp.geometry("500x500")
      resp.title("Resposta")
      
      Energia_label = ttk.Label(resp, text=f"Energia Fundamental: {np.format_float_scientific(E, precision = 2, exp_digits = 1)}", font=("Arial", 16))
      Energia_label.place(relx=0.5, rely=0.1, anchor=ttk.CENTER)

      Energia_Inicial_label = ttk.Label(resp, text=f"Energia Incial: {np.format_float_scientific(Eni, precision = 2, exp_digits = 1)} Ev", font=("Arial", 16))
      Energia_Inicial_label.place(relx=0.5, rely=0.3, anchor=ttk.CENTER)

      Energia_Final_label = ttk.Label(resp, text=f"Energia Final: {np.format_float_scientific(Enf, precision = 2, exp_digits = 1)} Ev", font=("Arial", 16))
      Energia_Final_label.place(relx=0.5, rely=0.5, anchor=ttk.CENTER)

      Energia_Foton_label = ttk.Label(resp, text=f"Energia do Fóton: {np.format_float_scientific(En, precision = 2, exp_digits = 1)} Ev", font=("Arial", 16))
      Energia_Foton_label.place(relx=0.5, rely=0.7, anchor=ttk.CENTER)

      Comprimento_label = ttk.Label(resp, text=f"Comprimento da Onda: {np.format_float_scientific(λ, precision = 2, exp_digits = 1)} m", font=("Arial", 16))
      Comprimento_label.place(relx=0.5, rely=0.9, anchor=ttk.CENTER)

    btn_calc = ttk.Button(jnl, text="Calcular", font=("Arial", 10), command=(calc3d))
    btn_calc.place(relx= 0.70, rely=0.4 , anchor=ttk.CENTER)
  
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
    selecionado_d['text'] = "Dimesão: 1D"

def def2d():
    global cx 
    cx = "2d"
    selecionado_d['text'] = "Dimesão: 2D"

def def3d():
    global cx 
    cx = "3d"
    selecionado_d['text'] = "Dimesão: 3D"


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
