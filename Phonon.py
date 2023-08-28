#!/usr/bin/env python
# coding: utf-8

# In[3]:


###########################################################################################################################
'''
                                Phonon frequency in a Fibonacci quasicrystal
                                        
                                            
                                            
This program generates a 1D quasicrystal using the Fibonacci generation sequence and calculates its phonon frequencies.
The phonon frequency is then plotted against phonon mode.

The user can choose the number of Fibonacci generations to create, the atom mass values, and a whether to print the atom
chain and mass occurence ratio, providing the masses are unique. The ratio is used to show that the generated crystal is 
a Fibonacci quasicrystal as it tends to the Golden ratio as the number of generations/atom chain length increases.
The user can also choose to generate the plots featuring increasing mass ratio and generation number, as featured in the
report in order to show the fractal behaviour of the frequency.

All the functions are put first, before the main section.

Note: Program significantly slows down for generations larger than about 15.


'''
###########################################################################################################################


##################################################### Imports #############################################################


import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import linalg
import math


################################################### Functions #############################################################


#Calcualtes analytical solution of phonon frequency given the atom chain
def analytical_sol(chain,mode):
    K = 1
    m = len(mode)
    sol = [0] * m
    for i in range(m):       
        sol[i] = 2*np.sqrt(K/chain[0])*np.sin(mode[i]*math.pi/(2*m))
    return(sol)


#Given initial masses and generation number, creates crystal atom chain
def crystal(gen,M_A,M_B): 
    
    #Variables
    sequence = [] 
    chain = []
    tempA = 'A'
    tempB = 'B'


    #Generate Fibonacci generations sequence as a string
    sequence.extend((tempA,tempB)) #Kick start generation with initial values
    for i in range(gen-2):
        nth = tempB + tempA
        sequence.append(nth)
        tempA = tempB
        tempB = nth
    #print("\n\nGeneration sequence: \n",sequence)


    #Separate generations sequence into atom chain
    for i in range(len(sequence)):
        chain += sequence[i]    
    
    
    #Convert mass placeholders into input floats
    for i in range(len(chain)):
        chain[i] = M_A if chain[i] == 'A' else M_B
    chain = np.array(chain,dtype=float)

    return(chain)


#Function to calculate phonon frequency eigenvalues given inital atom masses and number of generations, K value 
#and condition to print the atom chain and mass occurence ratio for proving quasicrystal
def Frequency(chain,K,print_atom_ratio=False): 

    #Generate matrix A for Eigenvalue equation:  A * Y = Omega^2 * Y
    m = len(chain)
    I = np.identity(m)
    A = np.zeros((m,m))
    for i in range(m):
        A[i][i] = -2*(chain[i]**(.5))*(-K)/chain[i]
        if i < len(A)-1:
            A[i][i+1] = (chain[i+1]**(.5))*(-K)/chain[i]
        if i > 0:
            A[i][i-1] = (chain[i-1]**(.5))*(-K)/chain[i]


    #Calculate frequency eigenvalues: Omega^2
    Eigenvalues, Eigenvectors = sp.linalg.eig(A)
    #print(Eigenvalues)


    #Generate mode numbers to plot against frequency
    mode = np.arange(0,m)
    

    #Prove Quasi-Crystal by calculating occurence of each mass in the chain. ratio for long chain tends to Golden ratio
    if print_atom_ratio == True:
        values, counts = np.unique(chain, return_counts=True)
        if len(values)>1: # Sanity check for 0 mass
            if counts[1] <=counts[0]:
                print("\nAtom chain:",chain)
                print("\n\nRatio of mass occurence for ",gen," generations =\n",counts[0]/counts[1])
            else:
                print("\nAtom chain:",chain)
                print("\n\nRatio of mass occurence for ",gen," generations = \n",counts[1]/counts[0])
    
    return(Eigenvalues,mode)


#Takes user input and returns if it's a positive integer
def get_pos_int():
    while True: 
        try:
            number = int(input("Input number of generations as a positive integer\n"))
            assert(number > 0)
            break
        except:
            print("Input must be a postive integer")
    return(number)


#Takes user input and returns if it's a positive float
def get_pos_float():
    while True:   
        try:
            number = float(input("Input mass of atom as a positive number\n"))
            assert(number > 0)
            break
        except:
            print("Input must be a postive number\n")
    return(number)


#Takes user input 'y' or 'n' and returns True or False
def get_bool(check):
    while True:
        if check == 'Plot':
            user = input("Generate plots of increasing mass and generations? Enter y/n\n")
        if check == 'Print':
            user = input("Print mass occurence ratio and atom chain?\n")
        if user == 'y':
            return(True)
            break
        elif user == 'n':
            return(False)
            break
        else:
            print("Please eneter either y/n\n")


#Plots phonon frequency vs mode number
def Plot(y,x):
    plt.figure(figsize=(10,10))
    plt.scatter(x,sorted(np.sqrt(y)))
    plt.title('Quasicrystal phonon frequency',fontsize=20)
    plt.ylabel('Phonon Frequency',fontsize=20)
    plt.xlabel('Mode',fontsize=20)
    plt.tick_params(axis = 'both',labelsize=18)
    plt.show()


#Plots the analytical soluton on top of the calculated solution 
def analytical_Plot(y1,x1,y2):
    plt.figure(figsize=(16,8))
    plt.scatter(x1,np.sqrt(sorted(y1)),label='Calculated')
    plt.plot(x1,y2,color='r', label='Analytical')
    plt.title('Phonon Frequency in a Quasicrystal',fontsize=20)
    plt.ylabel('Phonon Frequency',fontsize=20)
    plt.xlabel('Mode',fontsize=20)
    plt.tick_params(axis = 'both',labelsize=18)
    plt.legend(fontsize=20)
    plt.show()


#Generate increasing mass ratio and increasing generation plots
def ratio_Plot():
   
    for j in range(6,15): #Generations
        plt.figure(figsize=(16,8))
        for i in range(2,12): #Mass ratio
            f,g = Frequency(crystal(gen = j,M_A=1,M_B=i),K=1)
            plt.scatter(g,np.sqrt(sorted(f)),label = i)
        plt.title('Phonon frequency in a Quasi-Crystal',fontsize=20)
        plt.ylabel('Phonon Frequency',fontsize=20)
        plt.xlabel('Mode',fontsize=20)
        plt.tick_params(axis='both',labelsize=18)
        plt.legend(title='Mass ratio',fontsize=18)
        plt.show()


################################################# Main ####################################################################


#Get user input 
gen = get_pos_int()
mass2 = get_pos_float()
mass1 = get_pos_float()
plot_ratio = get_bool('Plot')


#Only print mass occurence ratio if masses are different
if mass1 != mass2:
    prnt = get_bool('Print')
else:
    prnt = False


#Calculate and return phonon frequencies; mode number; atom chain using user inputs
atom_chain = crystal(gen,mass1,mass2)
frequency, mode = Frequency(atom_chain,K=1,print_atom_ratio=prnt)


#Plot results, depending on inputs will produce analytical and caluclated solutions
if mass1 == 1 and mass2 ==1:
    sol = analytical_sol(atom_chain,mode)
    analytical_Plot(frequency,mode,sol)
else:
    Plot(frequency,mode)


#Mass ratio/fractal nature plots
if plot_ratio == True:
    ratio_Plot()

