#Należy aktywować poszczególne zadania w 'mainie', zmieniając wartości True/False

import numpy as np
import matplotlib . pyplot as plt
from scipy.integrate import solve_ivp
import math
import scipy.signal as sig

def step(t):
    if t>=0:
        return np.array([[1]])
    else:
        return np.array([[0]])

def noSignal(t):
    return np.array([[0]])

def sinSig2(t):
    return 0.1*np.sin(2*np.pi*2*t)

def sinSig0_65(t):
    return 0.1*np.sin(2*np.pi*0.65*t)

def sinSig0_2(t):
    return 0.1*np.sin(2*np.pi*0.2*t)

def sinSig10(t):
    return np.sin(10 * t)

def wahadlo(t, x, b, type):
    m=1
    g=9.81
    J=0.05
    l=1/2
    x = np.array([x]).T
    dx1=x[1,0]
    if type=='none':
        dx2=(noSignal(t)-m*g*l*np.sin(x[0,0])-b*x[1,0])/(m*l*l+J)
    elif type=='step':
        dx2 = (step(t) - m * g * l * np.sin(x[0, 0]) - b * x[1, 0]) / (m * l * l + J)
    elif type=='sin2':
        dx2 = (sinSig2(t) - m * g * l * np.sin(x[0, 0]) - b * x[1, 0]) / (m * l * l + J)
    elif type == 'sin0_65':
        dx2 = (sinSig0_65(t) - m * g * l * np.sin(x[0, 0]) - b * x[1, 0]) / (m * l * l + J)
    elif type == 'sin0_2':
        dx2 = (sinSig0_2(t) - m * g * l * np.sin(x[0, 0]) - b * x[1, 0]) / (m * l * l + J)
    dx=np.array([[dx1],[dx2]])
    return np.ndarray.tolist(dx.T[0])

def systemRLC(t,x,A,type):
    R=0.2
    L=0.1
    C=0.05
    x = np.array([x]).T
    i3=(0.25*x[1,0]/(5-x[1,0]))
    if type=='step':
        dx1=(A*step(t)-R*x[0,0]-x[1,0])/L
        dx2=(x[0,0]-i3)/C
    elif type=='sine':
        dx1 = (A * sinSig10(t) - R * x[0, 0] - x[1, 0]) / L
        dx2 = (x[0, 0] - i3) / C
    dx=np.array([[dx1],[dx2]])
    return np.ndarray.tolist(dx.T[0])

def zadanie2(active):
    if active:
        # inicjalizacja zmiennych
        m = 1
        k = 1
        b0 = 0
        b05 = 1 / 2
        b2 = 2
        # utworzenie macierzy
        A1 = np.array([[0, 1], [-k / m, -b0 / m]])
        A2 = np.array([[0, 1], [-k / m, -b05 / m]])
        A3 = np.array([[0, 1], [-k / m, -b2 / m]])
        B = np.array([[0], [1 / m]])
        C = np.array([1, 0])
        D = np.array([0])
        # przygotowanie sygnałów do symulacji
        t = np.linspace(0, 50, 5001)
        step_ = step(t)
        sine = np.sin(2 * t)
        # utworzenie systemów do symulacji
        system1 = sig.lti(A1, B, C, D)
        system2 = sig.lti(A2, B, C, D)
        system3 = sig.lti(A3, B, C, D)
        # symulacje
        ts1, ys1, xs1 = sig.lsim2(system1, step_, t)
        tsin1, ysin1, xsin1 = sig.lsim2(system1, sine, t)
        ts2, ys2, xs2 = sig.lsim2(system2, step_, t)
        tsin2, ysin2, xsin2 = sig.lsim2(system2, sine, t)
        ts3, ys3, xs3 = sig.lsim2(system3, step_, t)
        tsin3, ysin3, xsin3 = sig.lsim2(system3, sine, t)
        # wykresy
        plt.figure(1)
        plt.plot(t, ys1, label='Step system 1', color='r')
        plt.plot(t, ysin1, label='Sine system 1', color='k')
        plt.xlabel('Time')
        plt.ylabel('System 1')
        plt.legend()
        plt.figure(2)
        plt.plot(t, ys2, label='Step system 2', color='r')
        plt.plot(t, ysin2, label='Sine system 2', color='k')
        plt.xlabel('Time')
        plt.ylabel('System 1')
        plt.legend()
        plt.figure(3)
        plt.plot(t, ys3, label='Step system 3', color='r')
        plt.plot(t, ysin3, label='Sine system 3', color='k')
        plt.xlabel('Time')
        plt.ylabel('System 3')
        plt.legend()
        plt.show()
        # 2.3 - współczynnik tłumienia wpływa na oscylacje układu przy wymuszeniu skokiem jednostkowym,
        # tak jak sama nazwa wskazuje 'tłumik' i odpowiadający mu w równaniu 'współczynnik tłumienia
        # wpływa na to, jak szybko stłumione zostaną oscylacje w systemie,
        # im wyższy parametr b - tym szybciej stłumimy oscylacje
        # jeżeli b=0, oscylacje nie gasną

def zadanie3(active):
    if active:
        #3.1-3.3
        #inicjalizacja
        l = 1 / 2
        C=np.array([[l], [-l]])
        #solvery
        sol0 = solve_ivp(wahadlo, [0, 60], [np.pi / 2, 0], args=[0, 'none'], rtol=1e-10)
        x0=np.array(sol0.y[0,:])
        y0=C*x0
        sol1_10 = solve_ivp(wahadlo, [0, 60], [np.pi / 2, 0], args=[1 / 10, 'none'], rtol=1e-10)
        x1_10 = np.array(sol1_10.y[0, :])
        y1_10 = C * x1_10
        sol1_2 = solve_ivp(wahadlo, [0, 60], [np.pi / 2, 0], args=[1 / 2, 'none'], rtol=1e-10)
        x1_2 = np.array(sol1_2.y[0, :])
        y1_2 = C * x1_2
        #wykresy
        #
        plt.figure(0)
        plt.title('Distance - b=0')
        plt.plot(sol0.t,y0[0,:],label='X')
        plt.plot(sol0.t,y0[1,:],label='Y')
        plt.legend()
        #
        plt.figure(1)
        plt.title('Rodzina charakterystyk fazowych')
        plt.plot(sol0.y[0,:],sol0.y[1,:],label='b=0')
        plt.plot(sol1_10.y[0, :], sol1_10.y[1, :], label='b=1/10')
        plt.plot(sol1_2.y[0, :], sol1_2.y[1, :], label='b=1/2')
        plt.xlabel('Angle')
        plt.ylabel('w')
        plt.legend()
        #3.4
        #
        plt.figure(2)
        plt.title('Kąt')
        plt.plot(sol0.t,sol0.y[0,:],label='b=0')
        plt.plot(sol1_10.t, sol1_10.y[0, :], label='b=1/10')
        plt.plot(sol1_2.t, sol1_2.y[0, :], label='b=1/2')
        plt.xlabel('Czas')
        plt.ylabel('Kąt')
        #okres drgań na oko 1.81s
        #TODO 3.4 obliczenia?
        #3.5
        solSin2 = solve_ivp(wahadlo, [0, 60], [0, 0], args=[1/10, 'sin2'], rtol=1e-10)
        solSin0_65 = solve_ivp(wahadlo, [0, 60], [0, 0], args=[1 / 10, 'sin0_65'], rtol=1e-10)
        solSin0_2 = solve_ivp(wahadlo, [0, 60], [0, 0], args=[1 / 10, 'sin0_2'], rtol=1e-10)
        x2=np.array(solSin2.y[0,:])
        x0_65=np.array(solSin0_65.y[0, :])
        x0_2=np.array(solSin0_2.y[0, :])
        y2=C*x2
        y0_65 = C * x0_65
        y0_2 = C * x0_2
        plt.figure(3)
        plt.title('Wymuszenie sinusoidalne - l_X')
        plt.plot(solSin2.t,y2[0,:],label='f=2')
        plt.plot(solSin0_65.t, y0_65[0, :], label='f=0.65')
        plt.plot(solSin0_2.t, y0_2[0, :], label='f=0.2')
        plt.legend()
        plt.figure(4)
        plt.title('Wymuszenie sinusoidalne - l_Y')
        plt.plot(solSin2.t, y2[1, :], label='f=2')
        plt.plot(solSin0_65.t, y0_65[1, :], label='f=0.65')
        plt.plot(solSin0_2.t, y0_2[1, :], label='f=0.2')
        plt.legend()
        plt.show()
        #3.5 im mniejsza częstotliwość, tym szybszy okres danego wachadła,
        # a dodatkowo mniejsze odchylenie od odległości między środkiem masy a osią obrotu

def zadanie4(active):
    if active:
        #symulacje
        solStepN10 = solve_ivp(systemRLC, [0, 2], [0, 0], args=[-10, 'step'], rtol=1e-10)
        solStep2 = solve_ivp(systemRLC, [0, 2], [0, 0], args=[2, 'step'], rtol=1e-10)
        solStep5 = solve_ivp(systemRLC, [0, 2], [0, 0], args=[5, 'step'], rtol=1e-10)
        solStep10 = solve_ivp(systemRLC, [0, 2], [0, 0], args=[10, 'step'], rtol=1e-10)
        solSin2 = solve_ivp(systemRLC, [0, 2], [0, 0], args=[2, 'sine'], rtol=1e-10)
        solSin10 = solve_ivp(systemRLC, [0, 2], [0, 0], args=[10, 'sine'], rtol=1e-10)
        #wyciągnięcie wartości napięcia
        uStepN10=np.array(solStepN10.y[1,:])
        uStep2 = np.array(solStep2.y[1, :])
        uStep5 = np.array(solStep5.y[1, :])
        uStep10 = np.array(solStep10.y[1, :])
        uSin2 = np.array(solSin2.y[1, :])
        uSin10 = np.array(solSin10.y[1, :])
        #wyciągnięcie wartości prądu
        i1StepN10 = np.array(solStepN10.y[0, :])
        i1Step2 = np.array(solStep2.y[0, :])
        i1Step5 = np.array(solStep5.y[0, :])
        i1Step10 = np.array(solStep10.y[0, :])
        i1Sin2 = np.array(solSin2.y[0, :])
        i1Sin10 = np.array(solSin10.y[0, :])
        #wyliczenie wyjścia
        yStepN10 = i1StepN10 - ((0.25*uStepN10)/(5-uStepN10))
        yStep2 = i1Step2 - ((0.25*uStep2)/(5-uStep2))
        yStep5 = i1Step5 - ((0.25*uStep5)/(5-uStep5))
        yStep10 = i1Step10 - ((0.25*uStep10)/(5-uStep10))
        ySin2 = i1Sin2 - ((0.25*uSin2)/(5-uSin2))
        ySin10 = i1Sin10 - ((0.25*uSin10)/(5-uSin10))
        #wykresy
        plt.figure(0)
        plt.title('Prąd i_1 - step')
        plt.plot(solStepN10.t,i1StepN10,label='step, A=-10')
        plt.plot(solStep2.t, i1Step2, label='step, A=2')
        plt.plot(solStep5.t, i1Step5, label='step, A=5')
        plt.plot(solStep10.t, i1Step10, label='step, A=10')
        plt.legend()
        plt.figure(1)
        plt.title('Prąd i_1 - sin')
        plt.plot(solSin2.t, i1Sin2, label='sin, A=2')
        plt.plot(solSin10.t, i1Sin10, label='sin, A=10')
        plt.legend()
        #
        plt.figure(2)
        plt.title('Napięcie kondensatora - step')
        plt.plot(solStepN10.t, uStepN10, label='step, A=-10')
        plt.plot(solStep2.t, uStep2, label='step, A=2')
        plt.plot(solStep5.t, uStep5, label='step, A=5')
        plt.plot(solStep10.t, uStep10, label='step, A=10')
        plt.legend()
        plt.figure(3)
        plt.title('Napięcie kondensatora - sin')
        plt.plot(solSin2.t, uSin2, label='sin, A=2')
        plt.plot(solSin10.t, uSin10, label='sin, A=10')
        plt.legend()
        #
        plt.figure(4)
        plt.title('Prąd i_2 - step')
        plt.plot(solStepN10.t, yStepN10, label='step, A=-10')
        plt.plot(solStep2.t, yStep2, label='step, A=2')
        plt.plot(solStep5.t, yStep5, label='step, A=5')
        plt.plot(solStep10.t, yStep10, label='step, A=10')
        plt.legend()
        plt.figure(5)
        plt.title('Prąd i_2 - sin')
        plt.plot(solSin2.t, ySin2, label='sin, A=2')
        plt.plot(solSin10.t, ySin10, label='sin, A=10')
        plt.legend()
        #
        plt.show()

if __name__ == '__main__':
    zadanie2(False)
    zadanie3(False)
    zadanie4(True)