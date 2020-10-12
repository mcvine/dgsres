from __future__ import absolute_import
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
from .Ikeda_carpenter import Count as Count



def lamda(Ei):
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

        # b_cal = 1 / b0

        Ei_array=np.array([15,60,100,200,600])
        a_array=np.array([0.21850,0.35601,0.46473,0.71558,1.33170])
        b_array=np.array([0.03737,0.03931,0.04370,0.09717,0.19422])
        R_array=np.array([0.76500,0.73662,0.62021,0.31510,0.26689])
        t0_array = np.array([0.42578, 0.22595, 0.17405, 0.12609, 0.07622])

        A_fixed = 1
        s_fixed = 2


        a=QApplication([])

        mcvine = QFileDialog.getOpenFileNames()[0]

        # a_cal = interp1d(Ei_array, a_array, fill_value='interpolate', kind='cubic')(self.Ei)
        # b_cal = interp1d(Ei_array, b_array, fill_value='interpolate', kind='cubic')(self.Ei)
        # R_cal = interp1d(Ei_array, R_array, fill_value='interpolate', kind='cubic')(self.Ei)
        # t0_cal = interp1d(Ei_array, t0_array, fill_value='interpolate', kind='cubic')(self.Ei)
        if self.Ei == 30:
            a_cal=0.28175
            b_cal=0.03827
            R_cal=0.78498
            t0_cal=0.30687

        if self.Ei==130:
            a_cal = 0.528
            b_cal = 0.04910
            R_cal = 0.51
            t0_cal = 0.155

        if self.Ei==300:
            a_cal = 0.89
            b_cal = 0.15
            R_cal = 0.28
            t0_cal = 0.10
            # R_cal = 0.5  # 0.5 for Ei=30,

        # a_cal = 1. / (a0 + (self.lamda(self.Ei) * a1))

        no_Et_Files=len(mcvine)
        for i in range(no_Et_Files):
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

            # if self.Ei==30:
            #     if peak_position>20:
            #         s_fixed=8

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

            # t0_cal =1
            #
            # print (t_fixed)

            # if self.Ei==30:
            #     t_fixed = peak_position

            # R_cal=np.exp(-81.799/(k*(self.lamda(self.Ei)**2)))

            # if self.Ei == 130:
                # R_cal=0.5



            print(a_cal, b_cal, R_cal)

            # if self.Ei == 130:
            #         mod.set_param_hint('R', value=R_cal, min=0.0, max=1.0 )
            #         mod.set_param_hint('a', value=a_cal,  min=0.0)
            #         mod.set_param_hint('b', value=b_cal, min=0.0)
            #         mod.set_param_hint('s', value=s_fixed, min=0.0)
            #         mod.set_param_hint('t0', value=t_fixed, min=0.0)

            # if self.Ei == 30:
            #     mod.set_param_hint('R', value=R_cal, min=0.0, max=1.0)
            #     mod.set_param_hint('a', value=a_cal, min=0.0)
            #     mod.set_param_hint('b', value=b_cal, min=0.0)
            #     mod.set_param_hint('s', value=s_fixed, min=0.0)
            #     mod.set_param_hint('t0', value=t_fixed, min=0.0, vary=False)

            # if self.Ei == 30:
            #     mod.set_param_hint('R', value=R_cal, vary=False)
            #     mod.set_param_hint('a', value=a_cal, vary=False)
            #     mod.set_param_hint('b', value=b_cal, vary=False)
            #     mod.set_param_hint('s', value=s_fixed, min=0.0)
            #     mod.set_param_hint('t0', value=t_fixed, min=0.0)

            # if self.Ei == 300:
            if i ==0:
                mod.set_param_hint('R', value=R_cal, min=0.0, max=1.0 )
                mod.set_param_hint('a', value=a_cal)
                mod.set_param_hint('b', value=b_cal)
                mod.set_param_hint('s', value=s_fixed, min=0.0)
                mod.set_param_hint('t0', value=t0_cal, min=0.0)

            if i>0:
                mod.set_param_hint('R', value=R_cal, vary=False)
                mod.set_param_hint('a', value=a_cal,vary=False)
                mod.set_param_hint('b', value=b_cal,vary=False)
                mod.set_param_hint('s', value=s_fixed, min=0.0)
                mod.set_param_hint('t0', value=t0_cal, min=0.0, vary=False)

            mod.set_param_hint('p', value=p_val, min=0.0, vary=False)
            # a.exec_()
            pars = mod.make_params()


            result = mod.fit(RF_inpl, pars, tm=taxis, fit_kws={'nan_policy': 'omit'})

            RF_Fit.append(result.best_fit)

            residual = np.sqrt((np.sum((RF_inpl - result.best_fit) ** 2)) / len(EnerInt))

            error.append(residual)

            dely = mod.eval(pars, tm=taxis)

            R_cal0 = R_cal
            a_cal0 = a_cal
            b_cal0 = b_cal
            t0_cal0 = t0_cal

            R_cal=result.params['R'].value
            a_cal=result.params['a'].value
            b_cal=result.params['b'].value
            t0_cal=result.params['t0'].value

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

        np.savetxt('Ikeda_Carpenter Parameters for Ei= {}.xlsx'.format(self.Ei), np.c_[ peak_positionE, R_fit,a_fit, b_fit, s_fit, t0_fit ])

        return(RF_mcvine, RFE_mcvine,RF_intp,RFE_intp,RF_Fit, R_fit,  a_fit, b_fit, s_fit, t0_fit, p_fit, peak_positionE,error )

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

        for i, j in zip (Et_peak, range(len(Et_peak))):
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

def test ():
    model_30 = ModelFitting(30, 0.1, 3)
    param_30 = model_30.param()
    RF_mcvine_30, RFE_mcvine_30, RF_intp_30, RFE_intp_30, RF_Fit_30, R_30, a_30, b_30, s_30, t_30, p_30, peaks_30, error_30 = param_30

if __name__=='__main__':
    test()
