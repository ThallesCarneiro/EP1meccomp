import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

g = 9.8 #m/s^2
m_total = 1939 #kg
L = 2.95 #m
I = 1 #kgm^2
w = 10 #rad/s
mi = 0.42
beta = 0.02

m = 0.6*m_total
m1 = 0.4*m_total

#Relação entre as normais
F = beta*m*g
F1 = mi*m1*g

#Condições Iniciais

theta_d_i = 0
theta_i = 0.174533 # equivalente a 10graus
x_i = 0
x_d_i =0
condicoes_iniciais = np.array([x_i,x_d_i,theta_i,theta_d_i])

####################
#Implementação das Solucoes de EDOs
def edo(condicoes):
    #Recebe os valores de x1,x2,theta1 e theta2, retorna a array dy 
    x1 = condicoes[0]
    x2 = condicoes[1]
    theta1 = condicoes[2]
    theta2 = condicoes[3]

    theta2_d = (-2*I*w*(theta2**2)*m_total +m1*L*math.sin(theta1)*(F -F1*math.cos(theta1) +m1*L*(theta2**2)*math.cos(theta1)))/(m1*(L**2)*(m_total-m1*(math.sin(theta1)**2)))

    x2_d = (-F + F1*math.cos(theta1) -m1*L*(math.cos(theta1)*(theta2**2) +theta2_d*math.sin(theta1)))/(m_total)

    dy = np.array([x2,x2_d,theta2,theta2_d])

    return dy

def euler(c_i,t,dt):
    #Recebe as condições iniciais, tempo de integração e o passo de integração, retorna uma matriz 6*n de passos com a solução da EDO
    passos = math.floor(t/dt)
    solucao = np.array([[0.0]*4]*passos)
    aceleracoes =np.array([[0.0]*2]*passos)

    solucao[0] = c_i
    aceleracoes[0][0] = edo(c_i)[1]
    aceleracoes[0][1] = edo(c_i)[3]

    for n in range(1,passos):
        dy = edo(solucao[n-1])
        solucao[n] = solucao[n-1] + dt*dy
        aceleracoes[n][0] = dy[1]
        aceleracoes[n][1] = dy[3]

    return solucao, aceleracoes
       
def rk2(c_i,t,dt):
    passos = math.floor(t/dt)
    solucao = np.array([[0.0]*4]*passos)
    aceleracoes =np.array([[0.0]*2]*passos)

    solucao[0] = c_i
    aceleracoes[0][0] = edo(c_i)[1]
    aceleracoes[0][1] = edo(c_i)[3]    

    for n in range(1,passos):
        k1 = edo(solucao[n-1])
        k2 = edo(solucao[n-1] + (dt/2)*k1)

        solucao[n] = solucao[n-1] + dt*k2
        aceleracoes[n][0] = k2[1]
        aceleracoes[n][1] = k2[3]

    return solucao,aceleracoes

def rk4(c_i,t,dt):
    passos = math.floor(t/dt)
    solucao = np.array([[0.0]*4]*passos)
    aceleracoes =np.array([[0.0]*2]*passos)

    solucao[0] = c_i
    aceleracoes[0][0] = edo(c_i)[1]
    aceleracoes[0][1] = edo(c_i)[3]       

    for n in range(1,passos):
        k1 = edo(solucao[n-1])
        k2 = edo(solucao[n-1] + (dt/2)*k1)
        k3 = edo(solucao[n-1] + (dt/2)*k2)
        k4 = edo(solucao[n-1] + (dt)*k3)
        
        solucao[n] = solucao[n-1] + (dt/6)*(k1 + 2*k2 + 2*k3 + k4)
        aceleracoes[n][0] = k4[1]
        aceleracoes[n][1] = k4[3]

    return solucao,aceleracoes


####################
#Implementação da geração de gráficos

def graficos(solucao,aceleracoes,t,dt):
    fig, axis = plt.subplots(2,1)
    vel_lin,= axis[0].plot(np.arange(0.0,t,dt),solucao[:,1])
    ac_lin, =axis[0].plot(np.arange(0.0,t,dt),aceleracoes[:,0])
    axis[0].set(xlabel='tempo (s)',title=r' $\.x$ e $\ddotx$ com passo h =' + str(dt))
    axis[0].legend((vel_lin,ac_lin),(" $\.x$(m/s)"," $\ddotx$(m/s²)"))
    axis[0].grid()


    pos_ang,=axis[1].plot(np.arange(0.0,t,dt),solucao[:,2]*(180/math.pi))
    vel_ang = axis[1].plot(np.arange(0.0,t,dt),solucao[:,3]*(180/math.pi))
    axis[1].set(xlabel='tempo (s)',title=r' $\theta$ e $\.\theta$ com passo h =' + str(dt))
    axis[1].legend((vel_lin,ac_lin),('$\\theta$ (rad)','$\dot \\theta$ (rad/s)'))
    axis[1].grid()
    
    plt.show()


solucao,aceleracoes =rk2(condicoes_iniciais,20,0.01)
graficos(solucao,aceleracoes,20,0.01)