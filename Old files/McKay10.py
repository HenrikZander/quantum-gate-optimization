from qutip import *
import numpy as np
from functions import *
from plotting import *

delta0 = 10
H = 0.5 * delta0 * (sx1 * sx2 + sy1 * sy2)
initialState = tensor(excitedState,groundState,groundState)
timeStamps = np.linspace(0,1,500)

result = sesolve(H, initialState, timeStamps, [sz1, sz2]) 
plotExpect(result)
print("Process finished!")