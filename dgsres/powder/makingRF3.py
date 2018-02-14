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

def lamda(Ei): #
    """lambda calculation from initial energy

                Parameters
                ----------
                Ei :float
                    Incident Energy

                """
    Ej=Ei*1.6021*10**-22
    h_js=6.626*10**-34
    m_kg=1.674929*10**-27
    lam=h_js/np.sqrt(2*Ej*m_kg)
    return (lam)


def vFrmE(E):
    """velocity calculation from energy

                    Parameters
                    ----------
                    E :float
                        Energy

                    """
    Ej=E*1.6021*10**-22
    m=1.674929*10**-27
    v=np.sqrt((2.*Ej)/m)
    return(v)

def vsq_from_E(E):
    """square of the velocity calculation from energy

                        Parameters
                        ----------
                        E :float
                            Energy

                        """
    Ej=E*1.6021*10**-22
    m=1.674929*10**-27
    return (2.*Ej)/m


def t(l3,Ei,Et,Et_axis):
    """square of the velocity calculation from energy

                            Parameters
                            ----------
                            E :float
                                Energy

                            """
    Ef=Ei-Et
    T=(-(l3/vFrmE(Ef))+(l3/np.sqrt(vFrmE(Ei)**2-vsq_from_E(Et_axis))))*1e6
    return (T)


class ModelFitting():

    """making the Resolution Functions

            Parameters
            ----------
            Ei=float
                Incident energy

            Et_interpolation_space : float
                energy transfer spacing for the interpolation between 0 to Ei
            l3=float
                sample to detector distance

            """



    warnings.filterwarnings("ignore")



    def __init__(self, Ei, Et_interpolation_space,l3):

        self.Ei=Ei
        self.Et_interpolation_space=Et_interpolation_space
        self.l3=l3
        self.Count = Count
        self.lamda=lamda
        self.vFrmE=vFrmE
        self.t = t


    # plt.figure()
    def param(self):
        RF_mcvine=[]
        RFE_mcvine = []

        RF_intp=[]
        RFE_intp=[]
        RF_Fit=[]

        peak_positionE = []
        error=[]
        R_fit = []
        a_fit = []
        b_fit = []
        s_fit = []
        t0_fit = []
        A_fit = []
        p_fit = []

        a0 = 1.6
        a1 = 1.5
        b0 = 31.9
        k = 46.

        b_cal = 1 / b0


        A_fixed = 1
        s_fixed = 1


        a=QApplication([])

        mcvine = QFileDialog.getOpenFileNames()[0]

        no_Et_Files=len(mcvine)
        for i in xrange(no_Et_Files):
            mcvineSol=mcvine[i]

            print (mcvineSol)
            iqe = hh.load(mcvineSol)

            iqe.I[iqe.I != iqe.I] = 0

            ie = iqe.sum('Q')  # histogram of DOS

            ie.I[np.isnan(ie.I)] = 0

            ie.I = ie.I / np.sum(ie.I)

            RF_mcvine.append(ie.I)

            RFE_mcvine.append(ie.E)

            peak_position = ie.E[np.where(ie.I == max(ie.I))[0][0]]  # finding the position of the peak of Eenergy transfer thar are given from mcvine

            peak_positionE.append(peak_position)

            #### interpolation

            EnerInt = np.arange(min(ie.E), max(ie.E), self.Et_interpolation_space)

            RF_inpl = interp1d(ie.E, ie.I, fill_value='interpolate', kind='cubic')(EnerInt)

            RF_inpl[RF_inpl < 0] = 0

            RF_intp.append(RF_inpl)
            RFE_intp.append(EnerInt)

            ##model fitting
            taxis = self.t(self.l3, self.Ei, peak_position, EnerInt)

            mod = Model(self.Count)

            # RF_inpl = RF_inpl / np.sum(RF_inpl)



            p_val = np.sum(RF_inpl)
            #         RF_inter.append( RF_inpl(EnerInt ))

            t_fixed =1

            print (t_fixed)

            # if self.Ei==30:
            #     t_fixed = peak_position

            R_cal=np.exp(-81.799/(k*(self.lamda(self.Ei)**2)))

            if self.Ei == 130:
                R_cal=0.5

            # if self.Ei==30:
            #     R_cal = 0.5 # 0.5 for Ei=30,

            a_cal = 1. / (a0 + (self.lamda(self.Ei) * a1))

            if self.Ei == 130:
                    mod.set_param_hint('R', value=R_cal, min=0.0, max=1.0 )
                    mod.set_param_hint('a', value=a_cal,  min=0.0)
                    mod.set_param_hint('b', value=b_cal, min=0.0)
                    mod.set_param_hint('s', value=s_fixed, min=0.0)
                    mod.set_param_hint('t0', value=t_fixed, min=0.0)

            # if self.Ei == 30:
            #     mod.set_param_hint('R', value=R_cal, min=0.0, max=1.0)
            #     mod.set_param_hint('a', value=a_cal, min=0.0)
            #     mod.set_param_hint('b', value=b_cal, min=0.0)
            #     mod.set_param_hint('s', value=s_fixed, min=0.0)
            #     mod.set_param_hint('t0', value=t_fixed, min=0.0, vary=False)

            if self.Ei == 30:
                mod.set_param_hint('R', value=R_cal, vary=False)
                mod.set_param_hint('a', value=a_cal, vary=False)
                mod.set_param_hint('b', value=b_cal, vary=False)
                mod.set_param_hint('s', value=s_fixed, min=0.0)
                mod.set_param_hint('t0', value=t_fixed, min=0.0)

            if self.Ei == 300:
                    mod.set_param_hint('R', value=R_cal, vary=False )
                    mod.set_param_hint('a', value=a_cal,  vary=False)
                    mod.set_param_hint('b', value=b_cal, vary=False)
                    mod.set_param_hint('s', value=s_fixed, min=0.0)
                    mod.set_param_hint('t0', value=t_fixed, min=0.0)

            mod.set_param_hint('p', value=p_val, min=0.0, vary=False)
            # a.exec_()
            pars = mod.make_params()


            result = mod.fit(RF_inpl, pars, tm=taxis, fit_kws={'nan_policy': 'omit'})

            RF_Fit.append(result.best_fit)

            residual = np.sqrt((np.sum((RF_inpl - result.best_fit) ** 2)) / len(EnerInt))

            error.append(residual)

            dely = mod.eval(pars, tm=taxis)

            R_fit.append(result.params['R'].value)  # saving the parameters
            a_fit.append(result.params['a'].value)  # saving the parameters
            b_fit.append(result.params['b'].value)  # saving the parameters
            s_fit.append(result.params['s'].value)  # saving the parameters
            t0_fit.append(result.params['t0'].value)  # saving the parameters
            p_fit.append(result.params['p'].value)  # saving the parameters

        # print R_fit
        RF_mcvine =np.reshape(RF_mcvine, (no_Et_Files,-1))
        RFE_mcvine= np.reshape(RFE_mcvine, (no_Et_Files, -1))
        RF_intp= np.reshape(RF_intp, (no_Et_Files, -1))
        RFE_intp=np.reshape(RFE_intp, (no_Et_Files, -1))
        RF_Fit=np.reshape(RF_Fit, (no_Et_Files, -1))

        return(RF_mcvine , RFE_mcvine,RF_intp,RFE_intp,RF_Fit, R_fit, a_fit, b_fit, s_fit,t0_fit, p_fit, peak_positionE,error )

            # return ()

    def param_refining(self,param,peak_position):

        param_avg = np.average(param)
        del_param = np.array(np.where(param > param_avg))[0]
        param = np.delete(param, (del_param), axis=0)
        peak_param = np.delete(peak_position, (del_param), axis=0)
        return (param, peak_param)

    def param_interpol(self, param, peak, x):

        best_fit = np.poly1d(np.polyfit(peak, param, 1))(np.unique(peak))
        if best_fit[-1] < 0:
            best_fit = best_fit + np.abs(best_fit[-1]) + 0.1
        # best_fit=param
        try:
            param = interp1d(peak, best_fit, fill_value='extrapolate', kind='cubic')
        except ValueError:
            param = interp1d(peak, best_fit, fill_value='extrapolate', kind='linear')

        if param(x)<0.0:
            best_fit=best_fit+np.abs(param(x))+0.1

            try:
                param = interp1d(peak, best_fit, fill_value='extrapolate', kind='cubic')
            except ValueError:
                param = interp1d(peak, best_fit, fill_value='extrapolate', kind='linear')

        return (param)


    def makingRf (self, begin, end, space,R,a,b,s,t0):

        a0 = 1.6
        a1 = 1.5
        b0 = 31.9
        k = 46.

        b_cal = 1 / b0
        a_cal = 1. / (a0 + (self.lamda(self.Ei) * a1))
        R_cal = np.exp(-81.799 / (k * (self.lamda(self.Ei) ** 2)))


        Et_peak = np.arange(begin, end, space)

        Et_range_Haxis=np.arange(begin,end,space)


        time=np.zeros((len(Et_peak), len(Et_range_Haxis)))
        RF = np.zeros((len(Et_peak), len(Et_range_Haxis)))

        for i, j in zip (Et_peak, xrange(len(Et_peak))):
            tcheck = self.t(self.l3, self.Ei, i, Et_range_Haxis)
            time[j,:]=tcheck
            if self.Ei == 130:
                RF[j,:] = self.Count(tcheck, R(i), a(i), b(i), s(i), t0(i), 1)
            if self.Ei == 300:
                RF[j,:] = self.Count(tcheck, R_cal, a_cal, b_cal, s(i), t0(i), 1)
            if self.Ei == 30:
                RF[j,:] = self.Count(tcheck, R_cal, a_cal, b_cal, s(i), t0(i), 1)
        RF[RF < 0] = 0

        return(RF, time)

