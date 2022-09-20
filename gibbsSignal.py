#!/usr/bin/env python
# coding: utf-8

# In[ ]:


##task1
##import libraries numpy
from numpy import *
import matplotlib.pyplot as set1


# In[ ]:


#set size draft and font
set1.rcParams['figure.figsize']=[10,10]
set1.rcParams.update({'font.size':18})


# In[115]:


#set values and domain
dx = 0.001
L = pi
x = L*arange(-1+dx,1+dx,dx)
n = len(x)
nquart = int(floor(n/4))

#hat function 
f = zeros_like(x)
f[nquart:2*nquart]=(4/n)*arange(1,nquart+1)
f[2*nquart:3*nquart] = ones(nquart) - (4/n)*arange(0,nquart)

#draft figure line 
fig,fs1 = set1.subplots()
fs1.plot(x,f,'-',color='k',LineWidth=3)

#fourier series function
A0 = sum(f*ones_like(x))*dx
fFS = A0/2

A= zeros(20)
B= zeros(20) ## compute fourier series
for k in range(20):

    A[k]=sum(f*cos(pi*(k+1)*x/L))*dx
    B[k]=sum(f*sin(pi*(k+1)*x/L))*dx
    fFS = fFS + A[k]*cos((k+1)*pi*x/L)+ B[k]*sin((k+1)*pi*x/L)
    plt.title('Fourier Series approximate a continuous-time hat signal')
    fs1.plot(x,fFS,'-') ## show figure 


# In[ ]:





# In[101]:


fFS = (A0/2) * ones_like(f) #f*1 FFS amplitudes
maxmode = 100
F1 = zeros(maxmode)
F2 = zeros(maxmode) #variables
KL = zeros(maxmode)


# In[109]:


F1[0]=A0/2
KL[0]=linalg.norm(f-fFS)/linalg.norm(f) 

for k in range(1,maxmode): ##loop for mode 
    F1[k] = sum(f*cos(pi*k*x/L))*dx 
    F2[k] = sum(f*sin(pi*k*x/L))*dx 
    fFS = fFS + F1[m]*cos(k*pi*x/L)+ F2[k]*sin(k*pi*x/L) ##approximate cosines
    KL[k]=linalg.norm(f-fFS)/linalg.norm(f)
    
fig,axs = set1.subplots() #plot
axs.plot(arange(maxmode),F1,color='k',LineWidth=3) ## width
axs.semilogy(r,F1[r],'o',color='r',MarkerSize=7) ## points size


set1.show() ## plot show


# In[118]:


from numpy import * ## Task 2
import matplotlib.pyplot as ff1


# In[124]:



ff1.rcParams['figure.figsize']=[10,10] ## set size plot
ff1.rcParams.update({'font.size':15}) ## font size plot

dx = 0.01 ## variables domain
L = 2*pi
x = arange(0,L+dx,dx)
n = len(x)

nquart = int(floor(n/4))



# In[126]:


f = zeros_like(x) 
f[nquart:3*nquart]=1

A0 = sum(f*ones_like(x)*dx*2/L) 
fFS = A0/2 * ones_like(f)

for k in range(1,100): ## compute Gibbs phenomena
    A1 = sum(f*cos(2*pi*k*x/L))*dx *2/L
    B1 = sum(f*sin(2*pi*k*x/L))*dx*2/L
    fFS = fFS + A1*cos(2*k*pi*x/L) + B1*sin(2*k*pi*x/L)

ff1.plot(x,f,color='k',LineWidth=2) ## size of lines
ff1.plot(x,fFS,'-',color='r',LineWidth=1.2)
ff1.show() ## show plot


# In[ ]:




