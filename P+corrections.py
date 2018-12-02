import math
import numpy as np
import matplotlib.pyplot as plt
from math import gamma
from Changing_parameters import N,T,rho
from Tabular_values import a_0, Na, E_h
from Dirak_functions import integral_1_2, integral_3_2, integral_minus_1_2, igrek

from Atom_parameters import Atom_weight,z
from Cell import z_0, r_0, volume, theta
import scipy
from scipy import integrate
from State_functions import eta, rho_e

from Phi_to_V_transformation import mu

from working_progonka import PHI, X, RESULT3

from hi_function import HI, Diff_HI

from Pressure import P
from corrections import delta_P