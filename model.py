######################################################################################################################################################################

#      .oooooo.   oooo                  oooo                                                
#     d8P'  `Y8b  `888                  `888                                                
#    888           888 .oo.    .oooo.    888  ooo. .oo.  .oo.    .ooooo.  oooo d8b  .oooo.o 
#    888           888P"Y88b  `P  )88b   888  `888P"Y88bP"Y88b  d88' `88b `888""8P d88(  "8 
#    888           888   888   .oP"888   888   888   888   888  888ooo888  888     `"Y88b.  
#    `88b    ooo   888   888  d8(  888   888   888   888   888  888    .o  888     o.  )88b 
#     `Y8bood8P'  o888o o888o `Y888""8o o888o o888o o888o o888o `Y8bod8P' d888b    8""888P' 


# File name: model.py

# Author(s): Henrik Zander, Emil Ingelsten

# Date created: 29 April 2021

# Last modified: 29 April 2021

# Copyright 2021, Henrik Zander and Emil Ingelsten, All rights reserved.

######################################################################################################################################################################

from qutip import *
import numpy as np
from variables import *
from numba import njit

######################################################################################################################################################################
# Parameter function.

def getParameterBounds(maxAllowedGateTime=240):
    """
    This function gets the bounds for the different 
    parameters that the optimizer can change in the 
    simulation. The function assumes the associated parameters
    to have the shape: x = [Theta, delta, omegaPhi, omegaTB0, operationTime].
    ---------------------------------------------------------
    INPUT:
            maxAllowedGateTime (int) {Optional}: The maximum gate time that will be allowed as a solution from the optimizer.
    ---------------------------------------------------------
    OUTPUT:
            parameterBounds (array(tuples(int))): Array of tuples that each contain the associated upper and lower bounds for that parameter. 
    ---------------------------------------------------------
    """
    return [(-0.5,0.5),(0,0.25),(0,5),(27.5,47.5),(50,maxAllowedGateTime)]


######################################################################################################################################################################
# Functions that define the tunable bus flux signal.

@njit
def Phi(t, theta, delta, omegaphi):
    """
    This function calculates the magnetic flux when the AC part of the flux has 
    a constant amplitude. 
    """
    phi = theta + delta*np.cos(omegaphi*t)
    return phi


@njit
def tunableBus(t, theta, delta, omegaphi, omegatb0):
    """
    This function calculates the frequency for the tunable bus, for the case
    when the AC part of the flux has a constant amplitude.
    """
    oTB = omegatb0*np.sqrt(np.abs(np.cos(PI*Phi(t, theta, delta, omegaphi))))
    return oTB


def omegaTB(t, args):
    """
    Wrapper function for tunableBus that handles the variable assignments.
    """
    theta = args['theta']
    delta = args['delta']
    omegaphi = args['omegaphi']
    omegatb0 = args['omegatb0']
    omegaTBTh = args['omegaTBTh']
    return tunableBus(t, theta, delta, omegaphi, omegatb0) - omegaTBTh


@njit
def sinstep(x, x_min, x_max):
    """
    This function calculates the function value of a sinusoidal step function, 
    that starts at x_min and assumes the value 1 after x_max, at the point x.
    """
    x = (x - x_min)/(x_max - x_min)
    if x < 0:
        x = 0
    elif x > 1:
        x = 1
    result = 0.5 - 0.5*np.cos(np.pi*x)
    return result


def sinBox(t, operationTime):
    """
    This function calculates the value of a box function with a sinusoidal rise
    and fall, with a total length of operationTime, at the time t.
    A rise time of 25 ns is used.
    """
    tRise = 25
    tWait = operationTime - 2*tRise
    funVal = sinstep(t, 0, tRise) - sinstep(t, tWait + tRise, tWait + 2*tRise)
    return funVal


@njit
def PhiSinStep(t, theta, delta, omegaphi, sinBoxVal):
    """
    This function calculates the magnetic flux when the AC part of the flux has 
    a sinusoidal box envelope. 
    """
    phi = theta + sinBoxVal*delta*np.cos(omegaphi*t)
    return phi


@njit
def tunableBusSinStep(t, theta, delta, omegaphi, omegatb0, sinBoxVal):
    """
    This function calculates the frequency for the tunable bus, in the case
    where the AC part of the flux has a sinusoidal box envelope.
    """
    oTB = omegatb0*np.sqrt(np.abs(np.cos(PI*PhiSinStep(t, theta, delta, omegaphi, sinBoxVal))))
    return oTB


def omegaTBSinStep(t, args):
    """
    Wrapper function for tunableBusSinStep that handles the variable assignments.
    """
    theta = args['theta']
    delta = args['delta']
    omegaphi = args['omegaphi']
    omegatb0 = args['omegatb0']
    operationTime = args['operationTime']
    omegaTBTh = args['omegaTBTh']
    sinBoxVal = sinBox(t,operationTime)
    return tunableBusSinStep(t, theta, delta, omegaphi, omegatb0, sinBoxVal) - omegaTBTh


@njit
def coeffomegaTB(omegaTB0, Phi):
    """
    This function calculates the tunable bus frequency at a given flux Phi.
    """
    return omegaTB0*np.sqrt(np.abs(np.cos(np.pi*Phi)))


######################################################################################################################################################################
# Hamiltomnian function.


def getHamiltonian(x, N=2, eigEs=None, U_e=None, getBBHamiltonianComps=False, getEigenStatesBB = False, getEigenEnergies=False, sinStepHamiltonian=False):
    """
    This function creates the hamiltonian for the specified number 
    of energy levels. It also has the ability to return the hamiltonian
    in both the bare basis and the eigenbasis. The user can also choose 
    whether the function should return the hamiltonian using an AC flux
    with constant amplitude or with a sinusoidal box envelope.
    ---------------------------------------------------------
    INPUT:
            x (array(float)): An array containing the parameters needed to time evolve the hamiltonian.
            N (int) {Optional}: Specifies the number of energy levels that should be used in the hamiltonian. Defaults to 2 energy levels.
            eigEs (array(float)) {Optional}: An array of the eigenenergies that has been calculated for the bare basis hamiltonian.
            U_e (Qobj) {Optional}: A QuTiP quantum operator that changes the basis from the bare basis to the eigenbasis.
            sinStepHamiltonian (boolean) {Optional}: Chooses whether the AC part of the flux affecting the tunable coupler frequency should use a sinusoidal box envelope. Defaults to False.

        Set only ONE of these to True!:
            getBBHamiltonianComps (boolean) {Optional}: Get the hamiltonian in the bare basis without the time dependent coefficient. Defaults to False.
            getEigenStatesBB (boolean) {Optional}: Get the hamiltonian in the bare basis with the coefficent evaluated for the DC part of the flux. Defaults to False.
            getEigenEnergies (boolean) {Optional}: Get the hamiltonian in the bare basis as a function that in turn returns the hamiltonian for a specific DC flux. Defaults to False.
        
        If none of the above are True, eigEs and U_e has to be provided.
    ---------------------------------------------------------
    OUTPUT:
            Hamiltonian (array(Qobj)): The hamiltonian for the case specified by the input parameters.
    ---------------------------------------------------------
    """

    # Choose the number of energy levels that should be used in the creation of the hamiltonian components.
    if N==4:
        # Get the 4 energy level hamltonian components.
        H0BB = omegas[0]*ad4_1*a4_1 - (alphas[0]/2.0)*(1-ad4_1*a4_1)*ad4_1*a4_1 + omegas[1]*ad4_2*a4_2 - (alphas[1]/2.0)*(1-ad4_2*a4_2)*ad4_2*a4_2 - (alphas[2]/2.0)*(1-ad4_TB*a4_TB)*ad4_TB*a4_TB
        HiBB = gs[0]*(ad4_1 + a4_1)*(ad4_TB + a4_TB) + gs[1]*(ad4_2 + a4_2)*(ad4_TB + a4_TB)
        H1BB = ad4_TB*a4_TB
    elif N==3:
        # Get the 3 energy level hamltonian components.
        H0BB = omegas[0]*ad3_1*a3_1 - (alphas[0]/2.0)*(1-ad3_1*a3_1)*ad3_1*a3_1 + omegas[1]*ad3_2*a3_2 - (alphas[1]/2.0)*(1-ad3_2*a3_2)*ad3_2*a3_2 - (alphas[2]/2.0)*(1-ad3_TB*a3_TB)*ad3_TB*a3_TB
        HiBB = gs[0]*(ad3_1 + a3_1)*(ad3_TB + a3_TB) + gs[1]*(ad3_2 + a3_2)*(ad3_TB + a3_TB)
        H1BB = ad3_TB*a3_TB
    else:
        # Get the 2 energy level hamltonian components.
        H0BB = (-omegas[0]/2)*sz1 + (-omegas[1]/2)*sz2
        HiBB = gs[0]*(sp1*smTB + sm1*spTB) + gs[1]*(sp2*smTB + sm2*spTB)
        H1BB = (-1/2)*szTB
    
    # Choose which hamiltonian configuration is to be returned.
    if getBBHamiltonianComps:

        # Just return the bare basis hamiltonian components without the time dependent coefficient.
        return [H0BB,HiBB,H1BB]
    elif getEigenStatesBB:

        # Get the hamiltonian in the bare basis, where the time dependent coefficient is constant with a flux phi that is equal to the current theta parameter in x.
        # This hamiltonian is used to calculate the relevant eigenstates. 
        return H0BB + HiBB + coeffomegaTB(x[3],x[0])*H1BB
    elif getEigenEnergies:

        # Get the hamiltonian in the bare basis, that will be specified as a function so that the hamiltonian can be calculated for a dynamic constant flux phi.
        def hamiltonian(currentPhi):
            return H0BB + HiBB + coeffomegaTB(x[3],currentPhi)*H1BB
        return hamiltonian
    else:

        # The hamiltonian in the eigenbasis.
        HThEB = Qobj(np.diag(eigEs),dims=[[N,N,N],[N,N,N]])
        H1EB = U_e*H1BB*U_e.dag()

        if sinStepHamiltonian:
            # The time dependent coefficient will be specified with an AC part of the flux that will be sinusoidally modulated.
            return [HThEB, [H1EB, omegaTBSinStep]]
        else:
            # The time dependent coefficient will be specified with an AC part of the flux that will be constant.
            return [HThEB, [H1EB, omegaTB]]


######################################################################################################################################################################
# Helper functions for calculating the gate fidelity.


def getThetaEigenstates(x, H_const, H_omegaTB, omegaTBTh):
    """
    This function calculates the eigenstates and eigenenergies 
    for a hamiltonian in the bare basis, using the tunable
    bus frequency at a constant flux.
    """
    H = H_const + omegaTBTh * H_omegaTB
    return H.eigenstates()


def getEBUnitary(x,eigStsBB,nLevels,Dimension):
    """
    This function returns the unitary that is used for transforming
    from the bare basis into the eigenbasis.
    """
    # Construct U_e
    U_e = Qobj()
    for i in range(Dimension):
        U_e += Qobj(basis(Dimension,i),dims=[[nLevels,nLevels,nLevels],[1,1,1]]) * eigStsBB[1][i].dag()
    # NB: U_e is ordered based on eigenenergies
    return U_e


def getRFUnitary(Hrot,t):
    """
    This function returns the unitary that is used to transform
    a system into the rotating frame.
    """
    U_rf = (1j*Hrot*t).expm()
    return U_rf


######################################################################################################################################################################
# Gate fidelity function.

def getGateFidelity(x,N=2,wantiSWAP=False,wantCZ=False,wantI=False):
    """
    This function calculates the average gate fidelity for the 
    iSWAP and CZ quantum gates given a parameter set x.
    ---------------------------------------------------------
    INPUT:
            x (array(float)): An array containing the parameters needed to time evolve the hamiltonian.
            N (int) {Optional}: Specifies the amount of energy levels that should be used in the hamiltonian. Defaults to 2 energy levels.

        Only change ONE of these to True!:
            wantiSWAP (boolean) {Optional}: Specifies that the gate fidelity is to be calculated for the iSWAP gate. Defaults to False.
            wantCZ (boolean) {Optional}: Specifies that the gate fidelity is to be calculated for the CZ gate. Defaults to False.
            wantI (boolean) {Optional}: USED IN TESTING ONLY! Specifies that the gate fidelity is to be calculated for the identity gate. Defaults to False.
    ---------------------------------------------------------
    OUTPUT:
            F_avg (float): The average gate fidelity for the parameter set x and the choosen gate.  
    ---------------------------------------------------------
    """

    # Get all parts of the hamiltonian in the bare basis.
    HBBComps = getHamiltonian(x,N=N,getBBHamiltonianComps=True)

    # Given the number of considered energy levels for each qubit, the dimension of the combined tensor state is calculated.
    D = N**3
    # From here on, unless otherwise stated, every state is considered a tensor state.

    # We are especially interested in |000>, |010>, |100> and |110> since these make up the computational basis.
    if N==4:
        eigIndices = [0, 1, 2, 5]
    elif N==3:
        eigIndices = [0, 1, 2, 5]
    else:
        eigIndices = [0, 1, 2, 4]

    # Define a list r of eigenstates in the eigenbasis.
    r = []
    for ei in eigIndices:
        r.append(Qobj(basis(D,ei),dims=[[N,N,N],[1,1,1]]))

    # Calculate omegaTB at Phi = Theta
    omegaTBTh = coeffomegaTB(x[3],x[0])

    # Calculate eigenstates and eigenenergies in the bare basis at Phi = Theta
    eigStsBB = getThetaEigenstates(x, HBBComps[0]+HBBComps[1], HBBComps[2], omegaTBTh)

    # Get unitary for transformation into eigenbasis
    U_e = getEBUnitary(x, eigStsBB, N, D)

    # NB: r and U_e are ordered based on eigenenergies


    # Simulate evolution of eigenstates:

    # Determine the time stamps for which the evolution will be solved at.
    opTime = x[4]
    ts = np.linspace(0,opTime,int(3*opTime))
    # Calculate the eigenbasis hamiltonian
    HEB = getHamiltonian(x, N=N, eigEs=eigStsBB[0], U_e=U_e, sinStepHamiltonian=True)

    # Initialise a list c of the time-evolved eigenstates
    c = r
    # Calculate final states and store them in c
    for i in range(len(c)):
        output = sesolve(HEB, c[i], ts, args={'theta': x[0], 'delta': x[1], 'omegaphi': x[2], 'omegatb0': x[3], 'operationTime': x[4], 'omegaTBTh': omegaTBTh})
        c[i] = output.states[-1]
    # NB: The elements of each state in c is ordered based on eigenenergies

    # Calculate U_rf:
    U_rf = getRFUnitary(HEB[0] - U_e*HBBComps[1]*U_e.dag(), ts[-1])

    # Transform all the states in c to the rotating frame
    c_rf = U_rf * c

    # Calculate M-matrix such that M_ij = <r_i|c_j>_rf:
    # Initialize as a 4x4 zero nested list
    M = [[0,0,0,0], [0,0,0,0], [0,0,0,0], [0,0,0,0]]

    # Assign values to elements M_ij
    for i, ei in enumerate(eigIndices):
        for j, _ in enumerate(eigIndices):
            M[i][j] = c_rf[j][ei].item(0)

    if wantiSWAP:
        # Calculate phases (iSWAP):
        phi = np.angle(M[0][0])
        theta1 = np.angle(M[1][2]) + PI/2 - phi
        theta2 = np.angle(M[2][1]) + PI/2 - phi

        # Ideal iSWAP gate matrix (with phases):
        U = np.matrix([[np.exp(1j*phi), 0, 0, 0], [0, 0, np.exp(1j*(-PI/2 + theta1 + phi)), 0], [0, np.exp(1j*(-PI/2 + theta2 + phi)), 0, 0], [0, 0, 0, np.exp(1j*(theta1 + theta2 + phi))]])
    elif wantCZ:
        # Calculate phases (CZ):
        phi = np.angle(M[0][0])
        theta1 = np.angle(M[1][1]) - phi
        theta2 = np.angle(M[2][2]) - phi
        
        # Ideal CZ gate matrix (with phases):
        U = np.matrix([[np.exp(1j*phi), 0, 0, 0], [0, np.exp(1j*(theta1 + phi)), 0, 0], [0, 0, np.exp(1j*(theta2 + phi)), 0], [0, 0, 0, np.exp(1j*(PI + theta1 + theta2 + phi))]])
    elif wantI:
        # Calculate phases (I):
        phi = np.angle(M[0][0])
        theta1 = np.angle(M[1][1]) - phi
        theta2 = np.angle(M[2][2]) - phi
        
        # Ideal I gate matrix (with phases):
        U = np.matrix([[np.exp(1j*phi), 0, 0, 0], [0, np.exp(1j*(theta1 + phi)), 0, 0], [0, 0, np.exp(1j*(theta2 + phi)), 0], [0, 0, 0, np.exp(1j*(theta1 + theta2 + phi))]])

    # Change M's type to matrix to simplify calculating fidelity
    M = np.matrix(M)

    # print(M)

    # Calculate gate fidelity
    dimComputationalSubspace = 4
    F_avg = (np.absolute(np.trace(M*U.H))**2 + np.trace(M.H*M)) / (dimComputationalSubspace*(dimComputationalSubspace+1))

    return F_avg


######################################################################################################################################################################