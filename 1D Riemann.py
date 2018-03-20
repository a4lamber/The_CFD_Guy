"""
Description: 1D Riemann Problem numerical solution
             1D compressible, inviscid fluid
             it takes 2nd Order MacCormack Format


author:      Adam Zhang
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import xlwt
from tempfile import TemporaryFile
import sys

#properties
Gamma = 1.4 #diatomic
L = 2
TT = 0.4
Sf = 0.8
J = 1000

#global variables
#in order to avoid 0/0, a robust way of creation array is used
#kinda problematic, need to be fixed
U = np.ones((3,J+2),dtype=float)/1000000
Uf = np.ones((3,J+2),dtype=float)/1000000
Ef = np.ones((3,J+2),dtype=float)/1000000

#CFL condition
def CFL(U,dx):
    maxvel = 1e-100
    for i in range(1,J):
        u = U[1][i]/U[0][i]
        p = (Gamma-1)*(U[2][i]-0.5*U[0][i]*u*u)
        #print(p)
        vel = math.sqrt(abs(Gamma*p/U[0][i]))+ abs(u)
        if vel>maxvel:
            maxvel = vel
    return Sf*dx/maxvel

#Riemann IC
#classic sod shock tube conditions used here
rou1 = 1.0
u1 = 0.0
p1 = 1.0
rou2 = 0.125
u2 = 0.0
p2 = 0.1
dx =  L/J #spatial grid length

for i in range(0,int(J/2)):
    U[0][i] = rou1
    U[1][i] = rou1*u1
    U[2][i] = p1/(Gamma-1)+rou1*u1*u1/2

for j in range(int(J/2)+1,J+1):
    U[0][j] = rou2
    U[1][j] = rou2*u2
    U[2][j] = p2/(Gamma-1)+rou2*u2*u2/2

#Boundary conditions
def bound(U):
    #left boundary
    for k in range(0,3):
        U[k][0] = U[k][1]
    #right boundary
    for m in range(0,3):
        U[m][J+1] = U[m][J]

def U2E(U,E):
    u = U[1]/U[0]
    p = (Gamma-1)*(U[2]-0.5*U[1]*U[1]/U[0])

    E[0] = U[1]
    E[1] = U[0]*u*u +p
    E[2] = (U[2]+p)*u


def MacCormack_adam(U,Uf,Ef,dx,dt):
    r = dt/dx
    nu = 0.25
    for i in range(1,J):
        q = abs(abs(U[0][i+1]-U[0][i])-abs(U[0][i-1]-U[0][i])/(abs(U[0][i+1]-U[0][i])+abs(U[0][i-1]-U[0][i])+1e-100))
        for k in range(0,3):
            Ef[k][i] = U[k][i] +0.5*nu*q*(U[k][i+1]-2*U[k][i]+U[k][i-1])

    for k in range(0,3):
        for i in range(1,J):
            U[k][i] = Ef[k][i]

    for i in range(0,J+1):
        U2E(U[:,i],Ef[:,i])



    for i in range(0,J):
        for k in range(0,3):
            Uf[k][i] = U[k][i] - r*(Ef[k][i+1]-Ef[k][i])

    for i in range(0,J):
        U2E(Uf[:,i],Ef[:,i])

    for i in range(1,J):
        for k in range(0,3):
            U[k][i] = 0.5*(U[k][i]+Uf[k][i])-0.5*r*(Ef[k][i]-Ef[k][i-1])

sto_x = []
sto_rho = []
sto_u = []
sto_p = []
sto_e = []

def output2(U,dx):
    for i in range(0,J):
        rou = U[0][i]
        u = U[1][i]/rou
        p =(Gamma-1)*(U[2][i]-0.5*U[0][i]*u*u)
        e = U[2][i]

        sto_x.append(i*dx)
        sto_rho.append(rou)
        sto_u.append(u)
        sto_p.append(p)
        sto_e.append(e)



def output(U,dx):
    book = xlwt.Workbook()
    sheet1 = book.add_sheet("one")

    stuff = ["X","rho","u","p","e"]
    for j in range(0,5):
        sheet1.write(0,j,stuff[j])


    for i in range(0,J):
        rou = U[0][i]
        u = U[1][i]/rou
        p =(Gamma-1)*(U[2][i]-0.5*U[0][i]*u*u)
        e = U[2][i]
        sheet1.write(i+1,0,i*dx)
        sheet1.write(i+1,1,rou)
        sheet1.write(i+1,2,u)
        sheet1.write(i+1,3,p)
        sheet1.write(i+1,4,e)

    name = "Finally_god.xls"
    book.save(name)
    book.save(TemporaryFile())

#total time
T = 0
while (T < TT):
    dt = CFL(U,dx)
    T = T + dt
    MacCormack_adam(U,Uf,Ef,dx,dt)
    bound(U)
output2(U,dx)






plt.scatter(sto_x,sto_rho,label = "density distribution")
plt.scatter(sto_x,sto_u, label = "velocity distribution")
plt.scatter(sto_x,sto_p,label = "pressure distribution")
plt.legend()
plt.show()