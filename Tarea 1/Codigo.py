import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation


def Der(y,h):
    N=len(y);
    D=np.zeros(N);
    for i in range (0,N):
        if i==0:
            D[i]=(1/(2*h))*(-y[i+2]+4*y[i+1]-3*y[i]);
        if i==1:
            D[i]=(1/(2*h))*(y[i+1]-y[i-1]);
        if 1<i<N-2:
            D[i]=(1/(12*h))*(-y[i+2]+8*y[i+1]-8*y[i-1]+y[i-2]);
        if i==N-2:
            D[i]=(1/(2*h))*(y[i+1]-y[i-1]);
        if i== N-1:
            D[i]=(1/(2*h))*(3*y[i]-4*y[i-1]+y[i-2]);
    return D

def RK4(N,h,f,phi_0,x0):
    phi=np.zeros(N);
    x=np.zeros(N);
    for i in range(0, N):
        x[i]=x0+i*h;
        if (i == 0):
            phi[i]=phi_0;
        else:
            k1=f(x[i-1],phi[i-1]);
            k2=f(x[i-1]+h/2,phi[i-1]+k1*h/2);
            k3=f(x[i-1]+h/2,phi[i-1]+k2*h/2);
            k4=f(x[i-1]+h,phi[i-1]+k3*h);
            phi[i]=phi[i-1]+h*(k1+2*k2+2*k3+k4)/6;
    return phi,x

def RK4V(N,h,f,cond_0,t0):
    x=np.zeros([N,2]);
    t=np.zeros(N);
    for i in range(0, N):
        t[i]=t0+i*h;
        if (i == 0):
            for k in range(0,2):
                x[i,k]=cond_0[k];
        else:
            k1=f(t[i-1],x[i-1,:]);
            k2=f(t[i-1]+h/2,x[i-1,:]+k1*h/2);
            k3=f(t[i-1]+h/2,x[i-1,:]+k2*h/2);
            k4=f(t[i-1]+h,x[i-1,:]+k3*h);
            x[i,:]=x[i-1,:]+h*(k1+2*k2+2*k3+k4)/6;
    return x,t

def Dens_E1(phi,Dx_phi,Pot):
    N=len(phi);
    p=np.zeros(len(phi));
    for i in range(0,N):
        p[i]=(1/2)*(Dx_phi[i]**2)+Pot(phi[i]);  
    return p

def int(L,h):
    I=0;
    N=len(L)
    for i in range(0,N): 
        I+=L[i]*h;
    return I



# Ejercicio 4. Codigo para evolucionar en el tiempo la ecuacion de Klein-Gordon
def Dif_Fin_0(L,Nx,Nt,tf,U0,dU0,Pot,DerPot,CondicionF,Name1,Name2):
    dx=L/Nx;
    x0=-L/2;
    xf=L/2;
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
        
    #Condiciones de frotera
    #for j in range (2,Nt+1):
        #U[0,j]=CondicionF(x0,t[j]);
        #U[Nx,j]=CondicionF(xf,t[j]);
        
    u=np.zeros(Nx+1);    
    for i in range (1,Nx):
        u[i]=U[i,0]+dt*dU0(xi[i])+(1/2)*(C**2)*(U[i+1,0]-2*U[i,0]+U[i-1,0])+(1/2)*(dt**2)*(DerPot(U[i,0]));  
    u[0]=CondicionF(x0,t[1]);
    u[-1]=CondicionF(xf,t[1]);
    
    #u=np.zeros([Nx+1,Nt+1]);
    
    #Solucion por diferencias finitas

    for j in range (0,Nt):
        for i in range (1,Nx):
            if j==0:
                U[i,j+1]=-U[i,0]+2*u[i]+(C**2)*(u[i+1]-2*u[i]+u[i-1])+(dt**2)*(DerPot(u[i])); #-U[i,j-1]+2*U[i,j]+(C**2)*(U[i+1,j]-2*U[i,j]+U[i-1,j])+(dt**2)*(DerPot(U[i,j]));
            else:
                U[i,j+1]=-U[i,j-1]+2*U[i,j]+(C**2)*(U[i+1,j]-2*U[i,j]+U[i-1,j])+(dt**2)*(DerPot(U[i,j]));
        U[0,j+1]=CondicionF(x0,t[j]);
        U[-1,j+1]=CondicionF(xf,t[j]);
    
    #for i in range(1,Nx):
        #for j in range(0,Nt):
            #if j==0:
                #U[i,j+1]=U[i,j]+dt*dU0(xi[i])+(1/2)*(C**2)*(U[i+1,j]-2*U[i,j]+U[i-1,j])+(1/2)*(dt**2)*(DerPot(U[i,j])); 
            #else:
                #U[i,j+1]=-U[i,j-1]+2*U[i,j]+(C**2)*(U[i+1,j]-2*U[i,j]+U[i-1,j])+(dt**2)*(DerPot(U[i,j]));  

    p=np.zeros([Nx+1,Nt+1]);
    phi_x=np.zeros([Nx+1,Nt+1]);
    #phi_t=np.zeros([Nx+1,Nt+1]);
    
    #Derivada respecto a x:
    for j in range(0,Nt+1):
        phi_x[:,j]=Der(U[:,j],0.1);
        
    # Derivada respecto a t
    #for i in range(0,Nx+1):
        #phi_t[i,:]=Der(U[i,:],0.1);
        
    #Densidad de energia
    for i in range(0,Nx+1):
        for j in range(0,Nt+1):
            p[i,j]=(1/2)*(phi_x[i,j])**2+Pot(U[i,j]);
    
    #for i in range(0,Nx+1):
        #print(U[i,:],"Sol\n");
        #print(phi_t[i,:],"Der,t\n");
    
    #print(phi_x,"\n\n")
    #print(phi_t,"\n\n")
    #print(p,"\n\n")
    
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
    ax.set_ylim(-5,5)
    ax.set_title(Name1)
    plt.grid()
    ax.set_ylabel(r'$\phi(x,t)$')
    ax.set_xlabel(r'$x$')
    plt.show()
    ani.save('Solución.GIF', fps=200) 
    
    #Gráfica p(x,t)
    fig, ax = plt.subplots()
    line2 = ax.plot(xi, p[:,0],'bo--')[0]
    def update(frame):
        # for each frame, update the data stored on each artist.
        x = xi
        y = p[:,frame]
        # update the scatter plot:
        line2.set_xdata(x)
        line2.set_ydata(y)
        return line2
    ani = animation.FuncAnimation(fig=fig, func=update, frames = Nt)
    ax.set_xlim(xi[0],xi[-1])
    ax.set_ylim(0,15)
    ax.set_title(Name2)
    plt.grid()
    ax.set_ylabel(r'$ \rho(x,t)$')
    ax.set_xlabel(r'$x$')
    plt.show()
    ani.save('Densidad.gif', fps=200)
    return p
    
#E=np.zeros(10);
#for i in range(0,Nx+1):
    #for j in range(0,10):
        #E[j]+=p[i,j]*dx;
#print("\nEnergia= ",E)
def Energia(p,dx,Nt):
    E=np.zeros(Nt+1);
    for j in range (0,Nt+1):
        E[j]=int(p[:,j],dx);
    return E


# Codigo para resolver una EDP de dos variables U(t,x).
def Dif_Fin(L,Nx,Nt,tf,U0,dU0,CondicionF,Name):
    dx=L/Nx;
    x0=0;
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
        
    #Condiciones de frotera
    #for j in range (1,Nt+1):
        #U[0,j]=CondicionF(x0,t[j]);
        #U[Nx,j]=CondicionF(xf,t[j]);
        
    for i in range (1,Nx):
        U[i,1]=U[i,0]+dt*dU0(xi[i])+(1/2)*(C**2)*(U[i+1,0]-2*U[i,0]+U[i-1,0]);
    U[0,1]=CondicionF(x0,t[1]);
    U[-1,1]=CondicionF(xf,t[1]);
    
    #u=np.zeros([Nx+1,Nt+1]);
    
    for j in range (1,Nt):
        for i in range (1,Nx):
            U[i,j+1]=-U[i,j-1]+2*U[i,j]+(C**2)*(U[i+1,j]-2*U[i,j]+U[i-1,j]);
        U[0,j+1]=CondicionF(x0,t[j]);
        U[-1,j+1]=CondicionF(xf,t[j]);
    
    #Solucion por diferencias finitas
    #for i in range(1,Nx):
        #for j in range(0,Nt):
            #if j==0:
                #U[i,j+1]=U[i,j]+dt*dU0(xi[i])+(1/2)*(C**2)*(U[i+1,j]-2*U[i,j]+U[i-1,j])
            #else:
                #U[i,j+1]=-U[i,j-1]+2*U[i,j]+(C**2)*(U[i+1,j]-2*U[i,j]+U[i-1,j]) 
                
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
    ax.set_ylim(-2,2)
    ax.set_title(Name)
    plt.grid()
    ax.set_ylabel(r'$U(x,t)$')
    ax.set_xlabel(r'$x$')
    plt.show()
    ani.save('Solucion.GIF', fps=200)