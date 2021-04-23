from qutip import *
import numpy as np
from functions import *
from plotting import *
from variables import *
import qutip.logging_utils as logging
import qutip.control.pulseoptim as cpo
import datetime
from numba import njit

def getnLevels():
    return 2

def getD():
    return getnLevels()**3

@njit
def sinstep(x, x_min, x_max):
    x = (x - x_min)/(x_max - x_min)
    if x < 0:
        x = 0
    elif x > 1:
        x = 1
    result = 0.5 - 0.5*np.cos(np.pi*x)
    return result

@njit
def Phi(t, theta, delta, omegaphi):
    phi = theta + delta*np.cos(omegaphi*t)
    return phi


@njit
def tunableBus(t, theta, delta, omegaphi, omegatb0):
    oTB = omegatb0*np.sqrt(np.abs(np.cos(PI*Phi(t, theta, delta, omegaphi))))
    return oTB


def omegaTB(t, args):
    theta = args['theta']
    delta = args['delta']
    omegaphi = args['omegaphi']
    omegatb0 = args['omegatb0']
    return tunableBus(t, theta, delta, omegaphi, omegatb0)


def sinBox(t, operationTime):
    tRise = 25
    tWait = operationTime - 2*tRise
    funVal = sinstep(t, 0, tRise) - sinstep(t, tWait + tRise, tWait + 2*tRise)
    return funVal


@njit
def PhiSinStep(t, theta, delta, omegaphi, sinBoxVal):
    phi = theta + sinBoxVal*delta*np.cos(omegaphi*t)
    return phi


@njit
def tunableBusSinStep(t, theta, delta, omegaphi, omegatb0, sinBoxVal):
    oTB = omegatb0*np.sqrt(np.abs(np.cos(PI*PhiSinStep(t, theta, delta, omegaphi, sinBoxVal))))
    return oTB

 
def omegaTBSinStep(t, args):
    theta = args['theta']
    delta = args['delta']
    omegaphi = args['omegaphi']
    omegatb0 = args['omegatb0']
    operationTime = args['operationTime']
    sinBoxVal = sinBox(t,operationTime)
    return tunableBusSinStep(t, theta, delta, omegaphi, omegatb0, sinBoxVal)

@njit
def coeffomegaTB(omegaTB0, Phi):
    coeff = omegaTB0*np.sqrt(np.abs(np.cos(np.pi*Phi)))
    return coeff

def getProjectionOperator(eigenStateIndex=1):
    D = getD()
    n = getnLevels()
    pSt = Qobj(basis(D,eigenStateIndex),dims=[[n,n,n],[1,1,1]]) # 010 for eigenStateIndex = 1
    pOp = pSt * pSt.dag()
    return pOp


def getAllProjectionOperators():
    pOp1 = getProjectionOperator(2) # 100
    pOp2 = getProjectionOperator(1) # 010
    pOpTB = getProjectionOperator(3) # 001
    return [pOp1,pOp2,pOpTB]


def getHamiltonian(x, omegaTBTh=None, eigEs=None, U_e=None, getBBHamiltonianComps=False, getEigenStatesBB = False, getEigenEnergies=False, sinStepHamiltonian=False):
    #The format of x is the following: x = [Theta, delta, omegaPhi, omegaTB0, operationTime]
    H0BB = (-omegas[0]/2)*sz1 + (-omegas[1]/2)*sz2 + gs[0]*(sp1*smTB + sm1*spTB) + gs[1]*(sp2*smTB + sm2*spTB)
    H1BB = (-1/2)*szTB

    if getBBHamiltonianComps:
        return [H0BB,H1BB]
    elif getEigenStatesBB:
        return H0BB + coeffomegaTB(x[3],x[0])*H1BB
    elif getEigenEnergies:
        def hamiltonian(currentPhi):
            return H0BB + coeffomegaTB(x[3],currentPhi)*H1BB 
        return hamiltonian
    else:
        n = getnLevels()

        HThEB = Qobj(np.diag(eigEs),dims=[[n,n,n],[n,n,n]])
        H1EB = U_e*H1BB*U_e.dag()

        def omegaTBSinStep_new(t, args):
            return omegaTBSinStep(t, args) - omegaTBTh
        def omegaTB_new(t,args):
            return omegaTB(t,args) - omegaTBTh
        if sinStepHamiltonian:
            return [HThEB, [H1EB, omegaTBSinStep_new]]
        else:
            return [HThEB, [H1EB, omegaTB_new]]

'''
def getSStepHamiltonian(x,operationTime=300.0):
    #The format of x is the following: x = [Theta, deltamax, omegaPhi, omegaTB0]
    H0 = (-omegas[0]/2)*sz1 + (-omegas[1]/2)*sz2 + gs[0]*(sp1*smTB + sm1*spTB) + gs[1]*(sp2*smTB + sm2*spTB)
    H1 = (-1/2)*szTB
    
    tRise = 2.0
    tWait = operationTime - tRise
    smoothness = 1
    def smoothBox(t):
        return smoothstep(t, 0, tRise, smoothness) - smoothstep(t, tWait + tRise, tWait + 2*tRise, smoothness)
    def Phi(t):
        return x[0] + smoothBox(t) * x[1]*np.cos(x[2]*t)
    def omegaTB(t, args):
        return x[3]*np.sqrt(np.abs(np.cos(PI*Phi(t))))
    return [H0, [H1, omegaTB]]

# Denna är inpastead tillfälligt medan jag ser om jag kan få gate fidelity att funka:
def getSinStepHamiltonian(x,operationTime=300.0,tRise=25.0):
    #The format of x is the following: x = [Theta, deltamax, omegaPhi, omegaTB0]
    H0 = (-omegas[0]/2)*sz1 + (-omegas[1]/2)*sz2 + gs[0]*(sp1*smTB + sm1*spTB) + gs[1]*(sp2*smTB + sm2*spTB)
    H1 = (-1/2)*szTB
    
    tWait = operationTime - 2*tRise
    def sinBox(t):
        return sinstep(t, 0, tRise) - sinstep(t, tWait + tRise, tWait + 2*tRise)
    def Phi(t):
        return x[0] + sinBox(t) * x[1]*np.cos(x[2]*t)
    def omegaTB(t, args):
        return x[3]*np.sqrt(np.abs(np.cos(PI*Phi(t))))
    return [H0, [H1, omegaTB]]
'''

def getThetaEigenstates(x, H_const, H_omegaTB, omegaTBTh):
    H = H_const + omegaTBTh * H_omegaTB
    return H.eigenstates()


# Unitary for transforming from the bare basis into the eigenbasis:
def getEBUnitary(x,eigStsBB,nLevels):
    D = nLevels**3

    # Construct U_e
    U_e = Qobj()
    for i in range(D):
        U_e += Qobj(basis(D,i),dims=[[nLevels,nLevels,nLevels],[1,1,1]]) * eigStsBB[1][i].dag()
    # NB: U_e is ordered based on eigenenergies
    return U_e


# Unitary for transforming into the rotating frame
# NB: This is usable only when working in the eigenbasis
def getRFUnitary(HThEB,t):
    U_rf = (1j*HThEB*t).expm()
    return U_rf


def getGateFidelity(x,wantiSWAP=False,wantCZ=False,wantI=False):
    HBBComps = getHamiltonian(x,getBBHamiltonianComps=True)

    n = getnLevels()
    D = getD()

    # We are especially interested in |000>, |010>, |100> and |110>:
    eigIndices = [0, 1, 2, 4] 

    # Define a list r of eigenstates in the eigenbasis
    r = []
    for ei in eigIndices:
        r.append(Qobj(basis(D,ei),dims=[[n,n,n],[1,1,1]]))

    # Calculate omegaTB at Phi = Theta
    omegaTBTh = coeffomegaTB(x[3],x[0])

    # Calculate eigenstates and eigenenergies in the bare basis at Phi = Theta
    eigStsBB = getThetaEigenstates(x, HBBComps[0], HBBComps[1], omegaTBTh)

    # Get unitary for transformation into eigenbasis
    U_e = getEBUnitary(x, eigStsBB, n)

    # NB: r and U_e are ordered based on eigenenergies


    # Simulate evolution of eigenstates:
    # Calculate the eigenbasis hamiltonian
    opTime = x[4]
    ts = np.linspace(0,opTime,int(3*opTime))
    HEB = getHamiltonian(x, omegaTBTh=omegaTBTh, eigEs=eigStsBB[0], U_e=U_e, sinStepHamiltonian=True)

    # Initialise a list c of the time-evolved eigenstates
    c = r
    # Calculate final states and store them in c
    for i in range(len(c)):
        output = sesolve(HEB, c[i], ts, args={'theta': x[0], 'delta': x[1], 'omegaphi': x[2], 'omegatb0': x[3], 'operationTime': x[4]})
        c[i] = output.states[-1]
    # NB: c is ordered based on eigenenergies

    # Calculate U_rf:
    U_rf = getRFUnitary(HEB[0], ts[-1])

    # Transform c into the rotating frame
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
    N = 4
    F_avg = (np.absolute(np.trace(M*U.H))**2 + np.trace(M.H*M)) / (N*(N+1))

    return F_avg

def getEigenStates(x, hamiltonian):
    H = hamiltonian(x, getEigenStatesBB=True)
    eigenStates = H.eigenstates()
    return eigenStates[1][2], eigenStates[1][1] # 100, 010
    
def getInitialState(eigenStateIndex=2):
    D = getD()
    n = getnLevels()
    return Qobj(basis(D,eigenStateIndex),dims=[[n,n,n],[1,1,1]]) # 100 for eigenStateIndex = 2


def getInitialGuess():
    return [Theta, delta, omegaPhi, omegas[2], 100]


def getParameterBounds():
    #Format of x: x = [Theta, delta, omegaPhi, omegaTB0, operationTime]
    return [(-0.5,0.5),(0,0.25),(0,5),(27.5,60),(50,240)]


def timeEvolutionH1():
    H0 = (-omegas[0]/2)*sz1 + (-omegas[1]/2)*sz2 + gs[0]*(sp1*smTB + sm1*spTB) + gs[1]*(sp2*smTB + sm2*spTB)
    H1 = (-1/2)*szTB
    H = [H0, [H1, omegaTB]]
    timeStamps = np.linspace(0,3000,3000)
    initialState = tensor(excitedState,groundState,groundState) # The states of qubit 1, qubit 2 and tunable bus, in that order. (100)
    projectionOperators = [sz1, sz2, szTB] # Used to calculate different expected values for the state as a function of time. #THESE ARE NOT CORRECT CHANGE TO EMILS OPERATORS!!!!
    result = sesolve(H, initialState, timeStamps, projectionOperators)
    plotExpect(result)


"""def optimizePulseH1():
    logger = logging.get_logger()
    log_level = logging.WARN #logging.INFO
    example_name = 'McKay1PulseOptim'
    H0 = (-omegas[0]/2)*sz1 + (-omegas[1]/2)*sz2 + gs[0]*(sp1*smTB + sm1*spTB) + gs[1]*(sp2*smTB + sm2*spTB)
    H1 = [(-1/2)*szTB]
    n_ctrls = len(H1)
    initialState = tensor(excitedState,groundState,groundState)
    targetState = tensor(groundState,excitedState,groundState)
    n_ts = 500
    evo_time = 100
    infidelity_target = 1e-3
    max_iter = 1000
    max_wall_time = 200
    p_type = 'TRIANGLE'
    f_ext = None
    result = cpo.opt_pulse_crab_unitary(H0, H1, initialState, targetState, n_ts, evo_time, fid_err_targ=infidelity_target, max_iter=max_iter, max_wall_time=max_wall_time, init_coeff_scaling=5.0, num_coeffs=5, method_params={'xtol':1e-3}, guess_pulse_type=None, guess_pulse_action='modulate', out_file_ext=f_ext, log_level=log_level, gen_stats=True)
    result.stats.report()
    print("Final evolution\n{}\n".format(result.evo_full_final))
    print("********* Summary *****************")
    print("Final fidelity error {}".format(result.fid_err))
    print("Final gradient normal {}".format(result.grad_norm_final))
    print("Terminated due to {}".format(result.termination_reason))
    print("Number of iterations {}".format(result.num_iter))
    print("Completed in {} HH:MM:SS.US".format(datetime.timedelta(seconds=result.wall_time)))
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(2, 1, 1)
    ax1.set_title("Initial Control amps")
    ax1.set_ylabel("Control amplitude")
    ax1.step(result.time, np.hstack((result.initial_amps[:, 0], result.initial_amps[-1, 0])), where='post')
    ax2 = fig1.add_subplot(2, 1, 2)
    ax2.set_title("Optimised Control Amplitudes")
    ax2.set_xlabel("Time")
    ax2.set_ylabel("Control amplitude")
    ax2.step(result.time, np.hstack((result.final_amps[:, 0], result.final_amps[-1, 0])), where='post')
    plt.tight_layout()
    plt.show()
    
    
optimizePulseH1()"""

"""
def Phi(t):
    return Theta + delta*np.cos(omegaPhi*t)


def omegaTB(t, args):
    return omegas[2]*np.sqrt(np.abs(np.cos(PI*Phi(t))))
"""