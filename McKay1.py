from qutip import *
import numpy as np
from functions import *
from plotting import *
from variables import *


H0 = (-omegas[0]/2)*sz1 + (-omegas[1]/2)*sz2 + gs[0]*(sp1*smTB + sm1*spTB) + gs[1]*(sp2*smTB + sm2*spTB)
H1 = (-1/2)*szTB
H = [H0, [H1, omegaTB]]

times = np.linspace(0,500,500)

result = mesolve(H, )