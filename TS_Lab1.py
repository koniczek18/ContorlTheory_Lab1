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

def wahadloStep(t, x, b,type):
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
        sol0 = solve_ivp(wahadloStep, [0, 60], [np.pi / 2, 0], args=[0,'none'], rtol=1e-10)
        x0=np.array(sol0.y[0,:])
        y0=C*x0
        sol1_10 = solve_ivp(wahadloStep, [0, 60], [np.pi / 2, 0],args=[1/10,'none'], rtol=1e-10)
        x1_10 = np.array(sol1_10.y[0, :])
        y1_10 = C * x1_10
        sol1_2 = solve_ivp(wahadloStep, [0, 60], [np.pi / 2, 0],args=[1/2,'none'],rtol=1e-10)
        x1_2 = np.array(sol1_2.y[0, :])
        y1_2 = C * x1_2
        #wykresy
        plt.figure(0)
        plt.title('Distance - b=0')
        plt.plot(sol0.t,y0[0,:],label='X')
        plt.plot(sol0.t,y0[1,:],label='Y')
        plt.legend()
        plt.figure(1)
        plt.title('Rodzina charakterystyk fazowych')
        plt.plot(sol0.y[0,:],sol0.y[1,:],label='b=0')
        plt.plot(sol1_10.y[0, :], sol1_10.y[1, :], label='b=1/10')
        plt.plot(sol1_2.y[0, :], sol1_2.y[1, :], label='b=1/2')
        plt.xlabel('Angle')
        plt.ylabel('w')
        plt.legend()
        plt.show()


if __name__ == '__main__':
    zadanie2(False)
    zadanie3(True)
