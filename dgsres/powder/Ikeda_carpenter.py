from matplotlib import pyplot as plt
import numpy as np, histogram.hdf as hh, histogram as H
from scipy.interpolate import interp1d
from scipy.special import erf
from scipy.special import erfc
from scipy.optimize import minimize,curve_fit
from lmfit import Model
import sys
from PyQt5.QtWidgets import QApplication, QWidget, QInputDialog, QLineEdit, QFileDialog
from PyQt5.QtGui import QIcon
import warnings
import h5py


def Count(tm, R, a, b, s, t0, p):
    """Ikeda-Carpenter Model

                Parameters
                ----------
                tm=array_like
                    time axis
                a :float
                    parameter controlling the rise of the curve
                b :float
                    parameter controlling the fall of the curve
                s: floats
                    parameter controlling the widenning of the curve
                t0 : float
                    position of the curve
                p=float
                    Normalizing Factor

                """
    ul = ((1 / (np.sqrt(2) * s)) * (((11.6 * tm) / 16.6) - t0))

    vl_a = (ul + ((a / np.sqrt(2)) * ((s * 16.6) / (2. + 3.))))

    vl_b = (ul + ((b / np.sqrt(2)) * ((s * 16.6) / (2. + 3.))))

    C01l = (np.sqrt(np.pi / 2) * ((s * 16.6) / (2. + 3.)) * np.exp(vl_b ** 2 - ul ** 2) * erfc(vl_b))

    C02l = (np.sqrt(np.pi / 2) * ((s * 16.6) / (2. + 3.)) * np.exp(vl_a ** 2 - ul ** 2) * erfc(vl_a))

    C1l = (((s * 16.6) / (2. + 3.)) ** 2 * np.exp(vl_a ** 2 - ul ** 2) * (
    np.exp(-vl_a ** 2) - (np.sqrt(np.pi) * vl_a * erfc(vl_a))))

    C2l = (np.sqrt(2) * ((s * 16.6) / (2. + 3.)) ** 3 * np.exp(vl_a ** 2 - ul ** 2) * (
    (np.sqrt(np.pi) * (.5 + vl_a ** 2) * erfc(vl_a)) - (vl_a * np.exp(-vl_a ** 2))))

    C = ((1 / (16.6 * np.sqrt(2 * np.pi) * s)) * (((1 - R) * a ** 2 * C2l) + (
    (R * a ** 2 * b / (a - b) ** 3) * ((2 * C01l) - (((a - b) ** 2 * C2l) + (2 * (a - b) * C1l) + (2 * C02l))))))

    C[np.isnan(C)] = 0
    C[C < 0] = 0
    A = p / np.sum(C)
    return (A * C)

