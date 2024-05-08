#!/usr/bin/env python
# coding: utf-8

# In[20]:


#Import Section
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from IPython.display import HTML


# In[21]:


#Initializing Biological constants from Mourao 2011
#Note all chances have to be multiplied by dt, but dt=1
#Note as well um must be multiplied by the grid step size

#If we let one element be ~0.2um, and using the fact the paper used a cell of size 1000um^3 we need something a little bigger \
#than for the yeast cells. 50x50 would be sufficient

Conc=35#uM
griddim=10/50
DC=0.5*griddim*griddim*0.2 #um^2/s
scalefac=0.8
knuc1=0.1*griddim
knuc2=0.2*griddim

#Note paper I'm getting constants from assumed that the growth rate and shortening rate are dependent on the tubulin, 
#but this woud cause issues with my finite difference model. Plus they are both about the same at 0.192 and 0.218 respectively
kcat=10*10**(-6)*griddim
kcat1=0.06*griddim
kcat2=0.55*griddim
kres1=0.01*griddim
kres2=0.02*griddim


# In[22]:


#Make Mesh grid
L=1
pts=50 #Important for how big of an area we want
tot=1000 #Time for experiment 
dx=L/pts #Change in size for grid
dy=L/pts #Change in size for grid
dt=1 #Time steps
T=tot*dt #Total time
grid=np.zeros((pts,pts,T))
#Make DiffusionGrid
Dgrid=np.ones((L*pts,L*pts,T))


# In[23]:


#Making weird boundaries
#Want a 3/5 ratio of a to b to get close to our correct area
a=0.6*pts/2
b=0.8*pts/2
t=np.linspace(0,2*np.pi,100) #fine boundary for comparisons
xs=a*np.cos(t)
ys=b*np.sin(t)
#plt.plot(xs,ys)
#Making smaller area for the nucleation to occur in
n1=0.25*a
n2=0.25*a
t=np.linspace(0,2*np.pi,100)
n1x=n1*np.cos(t)
n2y=n2*np.sin(t)
xshift=0.4*a


# In[25]:


#Random Walk simulator
tempgrid=grid
MTs={}#Define dict which MTs live in
grows=0 #MT dict counter
Diffgrid=Dgrid
#Defining empty arrays which we add to later
#Note the MTs have a structure like a list of coordinates where the MT is, and then a state number, 
#all of which associated to a dict entry
x=0
y=0

avglengths=[]
numMt=[]
for i in range(0,T-1): #Walking through all time steps
    templengths=[]
    lengths=[]
    for m in range(0,pts):# Go through all grid points
        for n in range(0,pts):
            #Note here we shift the n points because our grid is written as a matrix not coordiantes
            if np.abs((m-pts/2)**2/n1**2)+np.abs((n-pts/2-xshift)**2/n2**2)<1 or np.abs((m-pts/2)**2/n1**2)+np.abs((n-pts/2+xshift)**2/n2**2)<1:
                if tempgrid[m,n,i]>=1:#If we are occupied, we skip the nucleation step
                    pass
                else:#Check for nucleation, add it to the growth counter
                    knucint=-(knuc1*Diffgrid[m,n,i]-knuc2)
                    if knucint>=np.random.random():
                        MTs["MT{0}".format(grows)]=np.array([[[m,n],1]])
                        grows=grows+1
    for q in range(len(MTs)):# Go through each MT, let it grow via random walk in biased direction
        mt=MTs["MT{0}".format(q)]#Pick individual MT
        if len(mt)==0: #If empty, skip. We do not delete dict entires
            pass
        else:
            if len(mt)==1::#Parsing dict entires
                x=mt[-1][0][0]
                y=mt[-1][0][1]
                state=mt[-1][1]
            else:    
                x=mt[-1][0][0]
                y=mt[-1][0][1]
                state=mt[-1][1]

            #Walking
            if state==1: #Note we want to restrict movement inside of cell, so if it would walk out, don't. 
                A=Diffgrid[x+1,y,i]
                if np.abs((x+1-pts/2)**2/a**2)+np.abs((y-pts/2)**2/b**2)>1:
                    A=0
                B=Diffgrid[x-1,y,i]
                if np.abs((x-1-pts/2)**2/a**2)+np.abs((y-pts/2)**2/b**2)>1:
                    B=0
                C=Diffgrid[x,y+1,i]
                if np.abs((x-pts/2)**2/a**2)+np.abs((y+1-pts/2)**2/b**2)>1:
                    C=0
                D=Diffgrid[x,y-1,i]
                if np.abs((x-pts/2)**2/a**2)+np.abs((y-1-pts/2)**2/b**2)>1:
                    D=0
                E=A+B+C+D
                trav=np.random.uniform(0,E)#Pick direction, add new coordinate to the list of all of the coordinatees of the MT
                if 0<=trav<A:
                    x=x+1
                    if np.abs((x-pts/2)**2/a**2)+np.abs((y-pts/2)**2/b**2)>1:
                        x=x-1
                    tempgrid[x,y,i+1:]=q #set grid point to a number denoting where the MT is
                    mt=np.vstack((mt,np.array([[x,y],1])))
                if A<=trav<B+A:
                    x=x-1
                    if np.abs((x-pts/2)**2/a**2)+np.abs((y-pts/2)**2/b**2)>1:
                        x=x+1
                    tempgrid[x,y,i+1:]=q
                    mt=np.vstack((mt,np.array([[x,y],1])))
                if B+A<=trav<C+B+A:
                    y=y+1
                    if np.abs((x-pts/2)**2/a**2)+np.abs((y-pts/2)**2/b**2)>1:
                        y=y-1
                    tempgrid[x,y,i+1:]=q
                    mt=np.vstack((mt,np.array([[x,y],1])))
                if C+B+A<=trav<D+C+B+A:
                    y=y-1
                    if np.abs((x-pts/2)**2/a**2)+np.abs((y-pts/2)**2/b**2)>1:
                        y=y+1
                    tempgrid[x,y,i+1:]=q
                    mt=np.vstack((mt,np.array([[x,y],1])))
                #kcatint=kcat
                kcatint=-Diffgrid[x,y,i]*kcat1+kcat2 #Setting variable catasrophe rate
                if kcatint>np.random.random():#Check for catasrophe after walking
                    if mt.ndim==2:
                        for j in range(len(mt)):
                            mt[j][1]=2#Set all of the states of the MTs to shrink (2)
                    else:
                        mt[1]=2
            if state==2: #if shirnking
                tempgrid[x,y,i+1:]=0#Set grid point to 0 since there is no MT
                mt=np.delete(mt,(-1),axis=0)#Delete last entry for MT coordinates
                Diffgrid[x,y,i+1]=Diffgrid[x,y,i]+Tubulinstart #Release tubulin to diffuse
                kresint=Diffgrid[x,y,i]*kres1+kres2#adjust rescue rate
                if kresint>np.random.random():#Check for rescue

                    if mt.ndim==2:
                        for j in range(len(mt)):
                            mt[j][1]=1#Set all of the states of the MTs to grow (1)
                    else:
                        mt[1]=1
        MTs["MT{0}".format(q)]=mt #Update Dict entry for each MT
    for m in range(0,pts-1): #Diffusion step
        for n in range(0,pts-1):
            if np.abs((m-pts/2)**2/a**2)+np.abs((n-pts/2)**2/b**2)>1: #Only want points inside of cell to diffuse
                grad=0
            else:
                grad=DC*((Diffgrid[m+1,n,i]-2*Diffgrid[m,n,i]+Diffgrid[m-1,n,i])/dx+(Diffgrid[m,n+1,i]-2*Diffgrid[m,n,i]+Diffgrid[m,n-1,i])/dy)
                if np.abs((m+1-pts/2)**2/a**2)+np.abs((n-pts/2)**2/b**2)>1:
                    grad=0   #Want to prevent boundary from acting like a source for diffusion
                if np.abs((m-1-pts/2)**2/a**2)+np.abs((n-pts/2)**2/b**2)>1:
                    grad=0
                if np.abs((m-pts/2)**2/a**2)+np.abs((n+1-pts/2)**2/b**2)>1:
                    grad=0
                if np.abs((m-pts/2)**2/a**2)+np.abs((n-1-pts/2)**2/b**2)>1:
                    grad=0
            Diffgrid[m,n,i+1]=Diffgrid[m,n,i]+grad*dt
    #Has to be done after this time step of the diffusion   
 #variables and subroutine to create the length and MT num graphs
    counter=0
    lens=[]
    for test in range(len(MTs)):
        mttemporary=MTs["MT{0}".format(test)]
        if len(mttemporary)==0:
            pass
        else:
            counter=counter+1
            lens.append(len(mttemporary))
    numMt.append(counter)
    avglengths.append(0.2*np.average(lens))
    
    
    for o in range(len(MTs)):#note we have to update the diffusion grid additions and subtractions after the diffusion step
#This is because we choose D[i+1]=D[i].
        mttemp=MTs["MT{0}".format(o)]
        if len(mttemp)==0:
            pass
        else:
            if len(mttemp)==1: #Parsing dict entires
                x=mttemp[-1][0][0]
                y=mttemp[-1][0][1]
                state=mttemp[-1][1]
            else:    
                x=mttemp[-1][0][0]
                y=mttemp[-1][0][1]
                state=mttemp[-1][1]
        if state==1:
            Diffgrid[x,y,i+1]=max(Diffgrid[x,y,i]-scalefac,0)
        else: #Releasing tubulin upon shrink
            Diffgrid[x,y,i+1]=min(1,Diffgrid[x,y,i]+scalefac)


# In[26]:


#Making Length Hist
plt.hist(avglengths)
plt.title('Histogram of MT Length')
plt.xlabel('Lengths (μm)')
#plt.savefig('Hist.png',dpi=600)


# In[27]:


#Making Number Graph
a=np.arange(0,len(numMt))
plt.plot(a,numMt)
plt.title('MT Number')
plt.xlabel('Time (s)')
plt.ylabel('Number of MTs')
#plt.savefig('Num.png',dpi=600)


# In[28]:


#Making Length Graph
a=np.arange(0,len(avglengths))
plt.plot(a,avglengths)
plt.title('Average MT Length Across Experiment')
plt.xlabel('Time (s)')
plt.ylabel('Avg Length of MTs (μm)')
#plt.savefig('Lengths.png',dpi=600)


# In[31]:


#Plotting and making a movie of the MTs over time
fig = plt.figure(figsize=(5,5))
#Note imshow thinks of matrices which start with (i,j) being row/column starting at the top left
#Regular coordinates count (i,j) being x,y or effectively column/row starting at the bottom right so we need to just how we visualize   
plt.plot(ys,xs)
plt.plot(n2y-xshift,n1x)
plt.plot(n2y+xshift,n1x)
#plt.plot(n1x,n2y-xshift)
#plt.plot(n1x,n2y+xshift)
im=plt.imshow(tempgrid[:,:,T-1],animated=True,extent=[-pts/2,pts/2,-pts/2,pts/2],cmap='magma')
plt.title('Simulation of Cell for T=%i Seconds' %T)
#fig.colorbar(im)
def animateframes(i):
    im.set_array(tempgrid[:,:,i])
    return [im]
anim=animation.FuncAnimation(fig,animateframes,frames=T,interval=100)
HTML(anim.to_jshtml())


# In[43]:


#Writing the movie to a gif
f = r"Microtubules.gif" 
writergif = animation.PillowWriter(fps=30) 
#anim.save(f, writer=writergif)


# In[32]:


#Plotting and making a movie of the MTs over time
fig2 = plt.figure(figsize=(5,5))
#plt.plot(xs,ys)
plt.plot(ys,xs)
plt.plot(n2y-xshift,n1x)
plt.plot(n2y+xshift,n1x)
im2=plt.imshow(Diffgrid[:,:,T-1],animated=True,extent=[-pts/2,pts/2,-pts/2,pts/2])
plt.title('Diffusion Grid for T=%i Seconds' %T)
def animateframes2(i):
    im2.set_array(Diffgrid[:,:,i])
    return [im2]
anim2=animation.FuncAnimation(fig2,animateframes2,frames=T,interval=100)
plt.show()
HTML(anim2.to_jshtml())


# In[ ]:


#Writing the movie to a gif
f = r"DiffusionGrid.gif" 
writergif = animation.PillowWriter(fps=30) 
#anim2.save(f, writer=writergif)


# In[ ]:




