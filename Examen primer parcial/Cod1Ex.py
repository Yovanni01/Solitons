# Examen 1er parcial. Problema 3b. Codigo para evolucionar en el tiempo un perfil inicial utilizando diferencias finitas

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation

# Codigo para resolver una EDP de dos variables U(t,x).
def Dif_Fin(L,Nx,Nt,tf,U0,f,Name,Name2):
    dx=2*L/Nx;
    x0=-L;
    xf=L;
    t0=0;
    dt=tf/Nt;
    C=dt/dx;
    U=np.zeros([Nx+1,Nt+1]);
    print("\nValor constante dt/dx= ",C);
    xi=np.zeros(Nx+1);
    t=np.zeros(Nt+1);
    
    # Construcción de la malla
    
    for i in range (0,Nx+1):
        xi[i]=x0+i*dx;
        for j in range(0,Nt+1):
            t[j]=t0+j*dt;
            
    #Condiciones iniciales
    
    for i in range (0,Nx+1):
        U[i,0]=U0(xi[i]);
    
    #Solucion por diferencias finitas
    
    for i in range (1,Nx):
        if i==1:
            U[i,1]=U[i,0]-(C/2)*(f(U[i+1,0])-f(U[i-1,0]))-(C/(2*dx**2))*(-3*U[i+4,0]+14*U[i+3,0]-24*U[i+2,0]+18*U[i+1,0]-5*U[i,0]);
        if 1<i<Nx-1:
            U[i,1]=U[i,0]-(C/2)*(f(U[i+1,0])-f(U[i-1,0]))-(C/(2*dx**2))*(U[i+2,0]-2*U[i+1,0]+2*U[i-1,0]-U[i-2,0]);
        if i==Nx-1:
            U[i,1]=U[i,0]-(C/2)*(f(U[i+1,0])-f(U[i-1,0]))-(C/(2*dx**2))*(3*U[i-4,0]-14*U[i-3,0]+24*U[i-2,0]-18*U[i-1,0]+5*U[i,0]);
    U[0,1]=0;
    U[-1,1]=0;
    
    #u=np.zeros([Nx+1,Nt+1]);
    
    for j in range (1,Nt):
        for i in range (1,Nx):
            if i==1:
                U[i,j+1]=U[i,j-1]-C*(f(U[i+1,j])-f(U[i-1,j]))-(C/(dx**2))*(-3*U[i+4,j]+14*U[i+3,j]-24*U[i+2,j]+18*U[i+1,j]-5*U[i,j]);
            if 1<i<Nx-1:
                U[i,j+1]=U[i,j-1]-C*(f(U[i+1,j])-f(U[i-1,j]))-(C/(dx**2))*(U[i+2,j]-2*U[i+1,j]+2*U[i-1,j]-U[i-2,j]);
            if i==Nx-1:
                U[i,j+1]=U[i,j-1]-C*(f(U[i+1,j])-f(U[i-1,j]))-(C/(dx**2))*(3*U[i-4,j]-14*U[i-3,j]+24*U[i-2,j]-18*U[i-1,j]+5*U[i,j]);
        U[0,j+1]=0;
        U[-1,j+1]=0;
        
    #Gráfica U(x)
    
    fig, ax = plt.subplots();
    ax.plot(xi,U[:,0],'b-');
    ax.set_title(Name2)
    #ax.set_xlim(0,7);
    #ax.set_ylim(0,2);
    ax.set_xlabel(r"$x$", fontsize=15, fontname='Times New Roman')
    ax.set_ylabel(r"$u(x)$", fontsize=15, fontname='Times New Roman')
    plt.grid()
    plt.show()
                
    #Gráfica U(x,t)
    
    fig, ax = plt.subplots()
    line2 = ax.plot(xi, U[:,0],'ro--')[0]
    def update(frame):
        # for each frame, update the data stored on each artist.
        x = xi
        y = U[:,frame]
        # update the scatter plot:
        line2.set_xdata(x)
        line2.set_ydata(y)
        return line2
    ani = animation.FuncAnimation(fig=fig, func=update, frames = Nt)
    ax.set_xlim(xi[0],xi[-1])
    #ax.set_ylim(-2,2)
    ax.set_title(Name)
    plt.grid()
    ax.set_ylabel(r'$u(x,t)$')
    ax.set_xlabel(r'$x$')
    plt.show()
    ani.save('Solucion.GIF', fps=200)
