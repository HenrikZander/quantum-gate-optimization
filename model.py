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

# Copyright 2021, Henrik Zander and Emil Ingelsten, All rights reserved.

######################################################################################################################################################################

from qutip import *
import numpy as np
from numba import njit

######################################################################################################################################################################
# Function for unpacking circuit parameters.


def unpackCircuitParameters(circuitData):
    frequencies = circuitData["frequencies"]
    anharmonicities = circuitData["anharmonicities"]
    couplings = circuitData["couplings"]

    return frequencies, anharmonicities, couplings

######################################################################################################################################################################
# Parameter function (DEPRECATED).
'''
def getParameterBounds(maxAllowedGateTime=240, wantTradCZ=False, wantTradiSWAP=False, wantCZ_20=False):
    """
    This function gets the bounds for the different
    parameters that the optimizer can change in the
    simulation. The function assumes the associated parameters
    to have the shape: x = [Theta, delta, omegaPhi, modulationTime].
    ---------------------------------------------------------
    INPUT:
            maxAllowedGateTime (int) {Optional}: The maximum gate time that
                will be allowed as a solution from the optimizer.
            wantTradCZ (boolean) {Optional}: If true, the parameter bounds
                are constricted so that |Theta| < phi_crossing and
                omegaPhi = omega_11<->20 ± 0.15 Grad/s (if wantCZ_20 is True), or
                omegaPhi = omega_11<->02 ± 0.15 Grad/s (if wantCZ_20 is False).
                Defaults to False.
            wantTradiSWAP (boolean) {Optional}: If true, the parameter bounds
                are constricted so that |Theta| < phi_crossing and
                so that omegaPhi is within 10% of omega_10<->01.
                Defaults to False.
            wantCZ_20 (boolean) {Optional}: If wantTradCZ, the value of this
                boolean determines which type of “traditional” CZ to search for
                (see wantTradCZ above). Defaults to False.
    ---------------------------------------------------------
    OUTPUT:
            parameterBounds (array(tuples(int))): Array of tuples that each
                contain the associated upper and lower bounds for that parameter.
    ---------------------------------------------------------
    """
    if (wantTradCZ):
        if (wantCZ_20):
            omegaPhiMin = omegaPhi_CZ_20 - 0.15
            omegaPhiMax = omegaPhi_CZ_20 + 0.15
        else:
            omegaPhiMin = omegaPhi_CZ_02 - 0.15
            omegaPhiMax = omegaPhi_CZ_02 + 0.15
        
        return [(-phi_crossing, phi_crossing), (0, 0.25), (omegaPhiMin, omegaPhiMax), (50, maxAllowedGateTime)]
    elif (wantTradiSWAP):

        omegaPhiMin = omegaPhi_iSWAP - 0.15
        omegaPhiMax = omegaPhi_iSWAP + 0.15
        
        return [(-phi_crossing, phi_crossing), (0, 0.25), (omegaPhiMin, omegaPhiMax), (50, maxAllowedGateTime)]
    else:
        return [(-0.5, 0.5), (0, 0.25), (0, 8), (50, maxAllowedGateTime)]
'''

######################################################################################################################################################################
# Functions used to create a sinBox-envelope.

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


def sinBox(t, operationTime, tRise):
    """
    This function calculates the value of a box function with a sinusoidal rise
    and fall, with a total length of operationTime, at the time t.
    A rise time of tRise (usually 25 ns) is used.
    """
    tWait = operationTime - 2*tRise
    funVal = sinstep(t, 0, tRise) - sinstep(t, tWait + tRise, tWait + 2*tRise)
    return funVal


######################################################################################################################################################################
# Functions that define the tunable bus flux signal.

@njit
def coeffomegaTB(omegaTB0, signalVal, useArccosSignal):
    """
    This function calculates the tunable bus frequency at a given flux Phi, 
        or at a given fraction of omegaTB0, when an arccos signal is used.
    """
    if useArccosSignal:
        return omegaTB0*signalVal
    else:
        return omegaTB0*np.sqrt(np.abs(np.cos(np.pi*signalVal)))

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Functions specific to cos signals

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
    oTB = omegatb0 * np.sqrt(np.abs(np.cos(np.pi*Phi(t, theta, delta, omegaphi))))
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
def PhiSinBox(t, theta, delta, omegaphi, sinBoxVal):
    """
    This function calculates the magnetic flux when the AC part of the flux has
    a sinusoidal box envelope.
    """
    phi = theta + sinBoxVal*delta*np.cos(omegaphi*t)
    return phi


@njit
def tunableBusSinBox(t, theta, delta, omegaphi, omegatb0, sinBoxVal):
    """
    This function calculates the frequency for the tunable bus, in the case
    where the AC part of the flux has a sinusoidal box envelope.
    """
    oTB = omegatb0 * np.sqrt(np.abs(np.cos(np.pi*PhiSinBox(t, theta, delta, omegaphi, sinBoxVal))))
    return oTB


def omegaTBSinBox(t, args):
    """
    Wrapper function for tunableBusSinBox that handles the variable assignments.
    """
    theta = args['theta']
    delta = args['delta']
    omegaphi = args['omegaphi']
    omegatb0 = args['omegatb0']
    operationTime = args['operationTime']
    omegaTBTh = args['omegaTBTh']
    tRise = args['riseTime']
    sinBoxVal = sinBox(t, operationTime, tRise)
    return tunableBusSinBox(t, theta, delta, omegaphi, omegatb0, sinBoxVal) - omegaTBTh


#---------------------------------------------------------------------------------------------------------------------------------------------------------------------
# Functions specific to arccos signals

@njit
def tunableBusArccosSinBox(t, dcAmplitude, acAmplitude, omegaphi, omegatb0, sinBoxVal):
    """
    This function calculates the frequency for the tunable bus, in the case
    where the AC part of the flux has a sinusoidal box envelope (and an arccos-signal is used).
    """
    oTB = omegatb0 * (dcAmplitude + sinBoxVal * acAmplitude * np.cos(omegaphi*t))
    return oTB


def omegaTBArccosSinBox(t, args):
    """
    Wrapper function for tunableBusArccosSinBox that handles the variable assignments.
    """
    dcAmplitude = args['dcAmplitude']
    acAmplitude = args['acAmplitude']
    omegaphi = args['omegaphi']
    omegatb0 = args['omegatb0']
    operationTime = args['operationTime']
    omegaTBTh = args['omegaTBTh']
    tRise = args['riseTime']
    sinBoxVal = sinBox(t, operationTime, tRise)
    return tunableBusArccosSinBox(t, dcAmplitude, acAmplitude, omegaphi, omegatb0, sinBoxVal) - omegaTBTh


######################################################################################################################################################################
# Hamiltonian function.


def getHamiltonian(x, N=2, eigEs=None, U_e=None, getBBHamiltonianComps=False, getEigenStatesBB=False, getEigenEnergies=False, sinBoxHamiltonian=False, useArccosSignal=False, circuitData=None):
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
            sinBoxHamiltonian (boolean) {Optional}: Chooses whether the AC part of the flux affecting the tunable coupler frequency should use a sinusoidal box envelope. Defaults to False.

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
    # Unpack the circuit data needed to fully specify the hamiltonian.
    omegas, alphas, gs = unpackCircuitParameters(circuitData)

    # Choose the number of energy levels that should be used in the creation of the hamiltonian components.
    if N > 2:
        # Definition of basic operators for an N-level system
        a = destroy(N)  # Annihilation operator for an N-level system
        ad = create(N)  # Creation operator for an N-level system
        I = qeye(N)  # Identity operator for an N-level system

        # The upgraded operators that are given by the tensor product between the different basic operators and the identity operator
        # Upgraded states are defined as qubit one, qubit two and tunable bus in that order
        a_1 = tensor(a, I, I)
        a_2 = tensor(I, a, I)
        a_TB = tensor(I, I, a)

        ad_1 = tensor(ad, I, I)
        ad_2 = tensor(I, ad, I)
        ad_TB = tensor(I, I, ad)

        # Get the N-energy level hamltonian components.
        H0BB = omegas[0]*ad_1*a_1 - (alphas[0]/2.0)*(1-ad_1*a_1)*ad_1*a_1 + omegas[1]*ad_2*a_2 - (alphas[1]/2.0)*(1-ad_2*a_2)*ad_2*a_2 - (alphas[2]/2.0)*(1-ad_TB*a_TB)*ad_TB*a_TB
        HiBB = gs[0]*(ad_1 + a_1)*(ad_TB + a_TB) + gs[1]*(ad_2 + a_2)*(ad_TB + a_TB)
        H1BB = ad_TB*a_TB
    else:
        # Definition of basic operators for a two-level system
        sp = sigmap()  # Raising operator for a two-level system
        sm = sigmam()  # Lowering operator for a two-level system
        sx = sigmax()  # Pauli sigma-x operator
        sy = sigmay()  # Pauli sigma-y operator
        sz = sigmaz()  # Pauli sigma-z operator
        I = qeye(2)  # Identity operator for a two-level system

        # The upgraded operators that are given by the tensor product between the different basic operators and the identity operator
        # Upgraded states are defined as qubit one, qubit two and tunable bus in that order
        sx1 = tensor(sx, I, I)
        sx2 = tensor(I, sx, I)
        sxTB = tensor(I, I, sx)

        sy1 = tensor(sy, I, I)
        sy2 = tensor(I, sy, I)
        syTB = tensor(I, I, sy)

        sz1 = tensor(sz, I, I)
        sz2 = tensor(I, sz, I)
        szTB = tensor(I, I, sz)

        sp1 = tensor(sp, I, I)
        sp2 = tensor(I, sp, I)
        spTB = tensor(I, I, sp)

        sm1 = tensor(sm, I, I)
        sm2 = tensor(I, sm, I)
        smTB = tensor(I, I, sm)

        # Get the 2 energy level hamltonian components.
        H0BB = (-omegas[0]/2)*sz1 + (-omegas[1]/2)*sz2
        HiBB = gs[0]*(sp1*smTB + sm1*spTB) + gs[1]*(sp2*smTB + sm2*spTB)
        H1BB = (-1/2)*szTB

    # Choose which hamiltonian configuration is to be returned.
    if getBBHamiltonianComps:

        # Just return the bare basis hamiltonian components without the time dependent coefficient.
        return [H0BB, HiBB, H1BB]
    elif getEigenStatesBB:

        # Get the hamiltonian in the bare basis, where the time dependent coefficient is constant with a flux phi that is equal to the current theta/DC-amplitude parameter in x.
        # This hamiltonian is used to calculate the relevant eigenstates.
        return H0BB + HiBB + coeffomegaTB(omegas[2], x[0], useArccosSignal=useArccosSignal) * H1BB
    elif getEigenEnergies:

        # Get the hamiltonian in the bare basis, specified as a function so that the hamiltonian can be calculated at a flux phi (or at a certain signal amplitude, when an arccos signal is used).
        def hamiltonian(signalVal):
            return H0BB + HiBB + coeffomegaTB(omegas[2], signalVal, useArccosSignal=useArccosSignal) * H1BB
        return hamiltonian
    else:

        # The hamiltonian in the eigenbasis.
        HThEB = Qobj(np.diag(eigEs), dims=[[N, N, N], [N, N, N]])
        H1EB = U_e*H1BB*U_e.dag()

        if sinBoxHamiltonian:
            # The AC part of the signal will be multiplied by a sinBox envelope.
            if useArccosSignal:
                omegaTBfun = omegaTBArccosSinBox
            else:
                omegaTBfun = omegaTBSinBox
        else:
            # The AC part of the signal will be multiplied by a constant envelope.
            omegaTBfun = omegaTB
        return [HThEB, [H1EB, omegaTBfun]]


######################################################################################################################################################################
# Helper functions for calculating the gate fidelity.


def getThetaEigenstates(H_const, H_omegaTB, omegaTBTh):
    """
    This function calculates the eigenstates and eigenenergies
    for a hamiltonian in the bare basis, using the tunable
    bus frequency at a constant flux.
    """
    H = H_const + omegaTBTh * H_omegaTB
    return H.eigenstates()


def getEBUnitary(eigStsBB, nLevels, Dimension):
    """
    This function returns the unitary that is used for transforming
    the bare basis hamiltonian into the eigenbasis hamiltonian.
    """
    # Construct U_e
    U_e = Qobj()
    for i in range(Dimension):
        U_e += Qobj(basis(Dimension, i), dims=[[nLevels, nLevels, nLevels], [1, 1, 1]]) * eigStsBB[1][i].dag()
    # NB: U_e is ordered based on eigenenergies
    return U_e


def getRFUnitary(Hrot, t):
    """
    This function returns the unitary that is used to transform
    a system into the rotating frame.
    """
    U_rf = (1j*Hrot*t).expm()
    return U_rf


# def getIndicesOld(N):
#     if N == 4:
#         eigIndices = [0, 1, 2, 5]
#     elif N == 3:
#         eigIndices = [0, 1, 2, 5]
#     else:
#         eigIndices = [0, 1, 2, 4]
#     return eigIndices


def getIndices(N, eigenstates, lowestOverlapAllowed=0.0):
    eigenIndices = []
    for q1 in range(2):
        for q2 in range(2):
            qTB = 0
            eigIndMaxOverlap = 0
            maxOverlap = 0
            for eigenstateIndex in range(len(eigenstates)):
                currentOverlap = Qobj(tensor(basis(N, q1), basis(N, q2), basis(N, qTB)), dims=[[N, N, N], [1, 1, 1]]).overlap(eigenstates[eigenstateIndex])
                if np.abs(currentOverlap) > np.abs(maxOverlap):
                    eigIndMaxOverlap = eigenstateIndex
                    maxOverlap = currentOverlap

            if (np.abs(maxOverlap) > lowestOverlapAllowed):
                eigenIndices.append(eigIndMaxOverlap)
            else:
                print(f'Bad Theta: At least one state in the computational subspace had a maximum eigenstate overlap below {lowestOverlapAllowed}')
                return None
    if len(eigenIndices) > len(set(eigenIndices)):
        print('Bad Theta: At least two states in the computational subspace got the same eigenindex')
        return None
    else:
        return eigenIndices


def eigenstateOrder(eigenstates, N):
    # dimension = N**3
    # assignedEigenstates = set()
    order = []
    for q1 in range(N):
        for q2 in range(N):
            for qTB in range(N):
                maxOverlap = (0, 0, 0)
                for eigenstateIndex in range(len(eigenstates)):
                    currentOverlap = Qobj(tensor(basis(N, q1), basis(N, q2), basis(N, qTB)), dims=[[N, N, N], [1, 1, 1]]).overlap(eigenstates[eigenstateIndex])
                    if np.abs(currentOverlap) > np.abs(maxOverlap[0]):
                        maxOverlap = (currentOverlap, (q1, q2, qTB), eigenstateIndex)

                #if (np.abs(maxOverlap[0]) > 0.9): # (not {maxOverlap[2]}.issubset(assignedEigenstates)) and
                    # assignedEigenstates.add(maxOverlap[2])
                #if len(eigenIndices) > len(set(eigenIndices)):
                order.append(maxOverlap)
    return order


def fidelityPostProcess(Hrot, c, ts, tIndices, eigIndices, iSWAP, SWAP, CZ, I):
    fidelities = []
    fidelityTimes = []

    for timeIndex in tIndices:
        currentStates = []
        for i in range(len(c)):
            currentStates.append(c[i][timeIndex])

        U_rf = getRFUnitary(Hrot, ts[timeIndex])

        # Transform all the states in c to the rotating frame
        c_rf = U_rf * currentStates

        # Calculate M-matrix such that M_ij = <r_i|c_j>_rf:
        # Initialize as a 4x4 zero nested list
        M = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]

        # Assign values to elements M_ij
        for i, ei in enumerate(eigIndices):
            for j, _ in enumerate(eigIndices):
                M[i][j] = c_rf[j][ei].item(0)

        if iSWAP:
            # Calculate phases (iSWAP):
            phi = np.angle(M[0][0])
            theta1 = np.angle(M[1][2]) + np.pi/2 - phi
            theta2 = np.angle(M[2][1]) + np.pi/2 - phi

            # Ideal iSWAP gate matrix (with phases):
            U = np.matrix([[np.exp(1j*phi), 0, 0, 0], [0, 0, np.exp(1j*(-np.pi/2 + theta1 + phi)), 0], [0, np.exp(1j*(-np.pi/2 + theta2 + phi)), 0, 0], [0, 0, 0, np.exp(1j*(theta1 + theta2 + phi))]])
        elif SWAP:
            # Calculate phases (SWAP):
            phi = np.angle(M[0][0])
            theta1 = np.angle(M[1][2]) - phi
            theta2 = np.angle(M[2][1]) - phi

            # Ideal SWAP gate matrix (with phases):
            U = np.matrix([[np.exp(1j*phi), 0, 0, 0], [0, 0, np.exp(1j*(theta1 + phi)), 0], [0, np.exp(1j*(theta2 + phi)), 0, 0], [0, 0, 0, np.exp(1j*(theta1 + theta2 + phi))]])
        elif CZ:
            # Calculate phases (CZ):
            phi = np.angle(M[0][0])
            theta1 = np.angle(M[1][1]) - phi
            theta2 = np.angle(M[2][2]) - phi

            # Ideal CZ gate matrix (with phases):
            U = np.matrix([[np.exp(1j*phi), 0, 0, 0], [0, np.exp(1j*(theta1 + phi)), 0, 0], [0, 0, np.exp(1j*(theta2 + phi)), 0], [0, 0, 0, np.exp(1j*(np.pi + theta1 + theta2 + phi))]])
        elif I:
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
        fidelities.append((np.absolute(np.trace(M*U.H))**2 + np.trace(M.H*M)).real/(dimComputationalSubspace*(dimComputationalSubspace+1)))
        fidelityTimes.append(ts[timeIndex])

    # print(fidelities)
    # print(fidelityTimes)
    return fidelities, fidelityTimes

######################################################################################################################################################################
# Gate fidelity function.


def getGateFidelity(x, N=2, iSWAP=False, SWAP=False, CZ=False, I=False, tIndices=[-76, -61, -23, -1], circuitData=None, riseTime=25.0, useArccosSignal=False, printResults=False):
    """
    This function calculates the gate fidelity for the
    iSWAP and CZ quantum gates given a parameter set x.
    ---------------------------------------------------------
    INPUT:
            x (array(float)): An array containing the parameters needed to time evolve the hamiltonian.
            N (int) {Optional}: Specifies the amount of energy levels that should be used in the hamiltonian. Defaults to 2 energy levels.
            tIndices (array(int)): Determines at which times to calculate the gate fidelity. Since the time evolution of the computational
                basis eigenstates is simulated for 75 extra time steps (each roughly 0.33 ns long) after the gate is applied, a time index
                of -76 corresponds to when the (oscillating part of) gate control signal is shut off, and a time index of -1 corresponds to
                the final simulated time, roughly 25 ns after the gate is applied. 

        Only change ONE of these to True!:
            iSWAP (boolean) {Optional}: Specifies that the gate fidelity is to be calculated for the iSWAP gate. Defaults to False.
            SWAP (boolean) {Optional}: Specifies that the gate fidelity is to be calculated for the SWAP gate. Defaults to False.
            CZ (boolean) {Optional}: Specifies that the gate fidelity is to be calculated for the CZ gate. Defaults to False.
            I (boolean) {Optional}: USED IN TESTING ONLY! Specifies that the gate fidelity is to be calculated for the identity gate.
                Defaults to False.
    ---------------------------------------------------------
    OUTPUT:
            fidelities, fidelityTimes (array(float), array(float)): The gate fidelity for the parameter set x and the selected gate,
                evaluated at times corresponding to the tIndices input, as well as those times.
    ---------------------------------------------------------
    """

    # Determine the time stamps for which the evolution will be solved.
    opTime = x[-1]
    nExtraSteps = 75

    ts1, stepSize = np.linspace(0, opTime, 3*int(opTime), retstep=True)
    ts2 = np.linspace(opTime + stepSize, opTime + nExtraSteps*stepSize, nExtraSteps)
    ts = np.append(ts1, ts2)

    # Signals with |A| + |B| > 1 are impossible in practice. Thus, we return bad fidelity at all times for such gates.
    if useArccosSignal and (abs(x[0]) + abs(x[1]) > 1):
        fidelities = []
        fidelityTimes = []
        for ti in tIndices:
            fidelities.append(1 - (abs(x[0]) + abs(x[1])))
            fidelityTimes.append(ts[ti])
            
        return fidelities, fidelityTimes

    
    # Change the units of omegaPhi from GHz to Grad/s
    omegaphi = x[2] * 2*np.pi
    
    # Unpack the circuit data needed to fully specify the hamiltonian.
    omegas, _, _ = unpackCircuitParameters(circuitData)

    # Get all parts of the hamiltonian in the bare basis.
    HBBComps = getHamiltonian(x, N=N, getBBHamiltonianComps=True, circuitData=circuitData, useArccosSignal=useArccosSignal)

    # Given the number of considered energy levels for each qubit, the dimension of the combined tensor state is calculated.
    D = N**3
    # From here on, unless otherwise stated, every state is considered a tensor state.

    # Calculate omegaTB at Phi = Theta
    # (or at the DC-amplitude x[0], in the case of an arccos signal)
    omegaTBTh = coeffomegaTB(omegas[2], x[0], useArccosSignal=useArccosSignal)

    # Calculate eigenstates and eigenenergies in the bare basis at Phi = Theta
    eigStsBB = getThetaEigenstates(HBBComps[0]+HBBComps[1], HBBComps[2], omegaTBTh)

    # We are especially interested in |000>, |010>, |100> and |110> since these make up the computational basis.
    # These states correspond (most closely) to the eigenstates with the following eigenindices:
    eigIndices = getIndices(N, eigStsBB[1])
    # print(eigIndices)

    if (eigIndices is not None):
        # Define a list r of eigenstates in the eigenbasis.
        r = []
        for ei in eigIndices:
            r.append(Qobj(basis(D, ei), dims=[[N, N, N], [1, 1, 1]]))

        # Get unitary for transformation into eigenbasis
        U_e = getEBUnitary(eigStsBB, N, D)

        # NB: r and U_e are ordered based on eigenenergies

    # Simulate evolution of eigenstates:

    # If eigenindices couldn't be generated, the function returns fidelity 0.2 at all examined timestamps.
    if eigIndices is None:
        fidelities = []
        fidelityTimes = []
        for ti in tIndices:
            fidelities.append(0.2)
            fidelityTimes.append(ts[ti])
        # print(fidelities)
        # print(fidelityTimes)
        return fidelities, fidelityTimes

    # Calculate the eigenbasis hamiltonian
    HEB = getHamiltonian(x, N=N, eigEs=eigStsBB[0], U_e=U_e, sinBoxHamiltonian=True, circuitData=circuitData, useArccosSignal=useArccosSignal)

    # Initialise a list c of the time-evolved eigenstates
    c = [stateToBeEvolved for stateToBeEvolved in r]
    
    # Calculate final states and store them in c
    if useArccosSignal:
        args = {'dcAmplitude': x[0], 'acAmplitude': x[1], 'omegaphi': omegaphi, 'omegatb0': omegas[2], 'operationTime': x[3], 'omegaTBTh': omegaTBTh, 'riseTime': riseTime}
    else:
        args = {'theta': x[0], 'delta': x[1], 'omegaphi': omegaphi, 'omegatb0': omegas[2], 'operationTime': x[3], 'omegaTBTh': omegaTBTh, 'riseTime': riseTime}
    
    for i in range(len(c)):
        output = sesolve(HEB, c[i], ts, args=args)
        c[i] = output.states
    # NB: The elements of each state in c is ordered based on eigenenergies

    # Pick out the rotating hamiltonian
    Hrot = np.array(HEB[0])
    Hrot[eigIndices[-1], eigIndices[-1]] = Hrot[eigIndices[-2], eigIndices[-2]] + Hrot[eigIndices[-3], eigIndices[-3]] - Hrot[eigIndices[0], eigIndices[0]]
    Hrot = Qobj(Hrot, dims=[[N, N, N], [N, N, N]])

    res = fidelityPostProcess(Hrot, c, ts, tIndices, eigIndices, iSWAP, SWAP, CZ, I)
    if printResults:
        print('--------------------------------------------------')
        #print(f'Energy levels: {N}')
        print(f'Solution: {x}')
        print(f'Fidelities: {res[0]}')
        print(f'Times [ns]: {res[1]}')
        #print('--------------------------------------------------')
        #print(res)
    return res


######################################################################################################################################################################
