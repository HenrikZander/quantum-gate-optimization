from qutip import *
import numpy as np
from functions import *
from plotting import *
from variables import *
import qutip.logging_utils as logging
import qutip.control.pulseoptim as cpo
import datetime
from numba import njit


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
def tunnableBus(t, theta, delta, omegaphi, omegatb0):
    oTB = omegatb0*np.sqrt(np.abs(np.cos(PI*Phi(t, theta, delta, omegaphi))))
    return oTB


def omegaTB(t, args):
    theta = args['theta']
    delta = args['delta']
    omegaphi = args['omegaphi']
    omegatb0 = args['omegatb0']
    return tunnableBus(t, theta, delta, omegaphi, omegatb0)


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
def tunnableBusSinStep(t, theta, delta, omegaphi, omegatb0, sinBoxVal):
    oTB = omegatb0*np.sqrt(np.abs(np.cos(PI*PhiSinStep(t, theta, delta, omegaphi, sinBoxVal))))
    return oTB

 
def omegaTBSinStep(t, args):
    theta = args['theta']
    delta = args['delta']
    omegaphi = args['omegaphi']
    omegatb0 = args['omegatb0']
    operationTime = args['operationTime']
    sinBoxVal = sinBox(t,operationTime)
    return tunnableBusSinStep(t, theta, delta, omegaphi, omegatb0, sinBoxVal)


def getProjectionOperators():
    pSt1 = tensor(eSt,gSt,gSt) # 100
    pOp1 = pSt1 * pSt1.dag()
    pSt2 = tensor(gSt,eSt,gSt) # 010
    pOp2 = pSt2 * pSt2.dag()
    return [pOp2]


def getAllProjectionOperators():
    pSt1 = tensor(eSt,gSt,gSt) # 100
    pOp1 = pSt1 * pSt1.dag()
    pSt2 = tensor(gSt,eSt,gSt) # 010
    pOp2 = pSt2 * pSt2.dag()
    pStTB = tensor(gSt,gSt,eSt) # 001
    pOpTB = pStTB * pStTB.dag()
    return [pOp1,pOp2,pOpTB]


def getHamiltonian(x, getEigenStates = False, getEigenEnergies=False, sinStepHamiltonian=False):
    #The format of x is the following: x = [Theta, delta, omegaPhi, omegaTB0]
    H0 = (-omegas[0]/2)*sz1 + (-omegas[1]/2)*sz2 + gs[0]*(sp1*smTB + sm1*spTB) + gs[1]*(sp2*smTB + sm2*spTB)
    H1 = (-1/2)*szTB
    if getEigenStates:
        return H0 + x[3]*np.sqrt(np.abs(np.cos(PI*x[0])))*H1
    elif getEigenEnergies:
        def hamiltonian(currentPhi):
            return H0 + x[3]*np.sqrt(np.abs(np.cos(PI*currentPhi)))*H1
        return hamiltonian
    elif sinStepHamiltonian:
        return [H0, [H1, omegaTBSinStep]]
    else:
        return [H0, [H1, omegaTB]]


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


def getEigenStates(x, hamiltonian):
    H = hamiltonian(x, getEigenStates=True)
    eigenStates = H.eigenstates()
    return eigenStates[1][2], eigenStates[1][1]
    
    
def getInitialState():
    return tensor(excitedState,groundState,groundState) #|100> ket


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