# -*- Python -*-
#
# Jiao Lin <jiao.lin@gmail.com>
#

# Ikeda Carpenter function

import numpy as np
from scipy.special import erfc

def resolution(E, Ei, E0, a, b, R, sigma=5., t0=0., geom=None):
    """compute reoslution function given energy transfer axis (E)
    and other parameters.
    The resolution function calculated should be centered around E0,
    the nominal energy transfer.

    E: energy transfer axis
    Ei: incident energy
    E0: nominal energy transfer
    a,b,R Ikeda Carpenter function parameters
    sigma: additional gaussian broadening
    t0: emission time
    geom: instrument geometry. Geom instance
    """
    from mcni.utils.conversion import e2v, SE2V
    source = Source(); source.a = a; source.b = b; source.R = R
    icg = ICG(source, sigma, geom, t0)
    Ef = Ei - E0
    vf = e2v(Ef); vi = e2v(Ei)
    l3 = geom.l3
    t = -l3/vf + l3/np.sqrt(vi**2-SE2V**2*E)
    t *= 1e6
    r = icg.ICG(t)
    r[r!=r] = 0 # t may overflow, for them r should be 0
    return r


class ICG:
    """
    class for an instance of an ikeda carpenter function convoluted with a gaussian
    """
    def __init__(self, source, sigma, geom, t0):
        self.source = source
        self.sigma = sigma
        self.geom = geom
        self.t0 = t0
    def ICG(self, t):
        l, l2, l3 = self.geom.l, self.geom.l2, self.geom.l3
        a, b, R = self.source.a, self.source.b, self.source.R
        sigma = self.sigma
        vmin_a = self.vmin(a,t); vmin_b = self.vmin(b, t)
        umin = self.umin(t)
        sqp = np.sqrt(np.pi); sq2 = np.sqrt(2)
        T1 = sigma*l/(l2+l3);  T2_a = np.exp(vmin_a**2 - umin**2)
        C0_a = sqp/sq2*T1*T2_a*erfc(vmin_a)
        C1_a = T1**2 * T2_a * (np.exp(-vmin_a*vmin_a) - sqp*vmin_a * erfc(vmin_a))
        C2_a = sq2 * T1**3 * T2_a * (sqp*(1./2+vmin_a**2)*erfc(vmin_a)-vmin_a*np.exp(-vmin_a**2))
        T2_b = np.exp(vmin_b**2 - umin**2)
        C0_b = sqp/sq2*T1*T2_b*erfc(vmin_b)
        return 1./l/sq2/sqp/sigma*((1-R)*a*a*C2_a + R*a*a*b/(a-b)**3*(2*C0_b-((a-b)**2*C2_a+2*(a-b)*C1_a+2*C0_a)))
    def umin(self, t):
        sigma = self.sigma
        geom = self.geom
        return 1./np.sqrt(2.)/sigma * (geom.l1/geom.l * t - self.t0)
    def vmin(self, x, t):
        sigma = self.sigma
        geom = self.geom
        return self.umin(t) + x/np.sqrt(2) * (sigma*geom.l/(geom.l2+geom.l3))

class Source:
    """
    """
    def __init__ (self,a=None,b=None,R=None):
        self.a = a
        self.b = b
        self.R = R

class Geom:
    """
    l1: moderator to fermi chopper
    l2: fermi chopper to sample
    l3: sample to detector pixel
    """
    def __init__(self, l1, l2, l3):
        self.l1 = l1
        self.l2 = l2
        self.l3 = l3
        self.l = l1 + l2 + l3
        return


def test1():
    source = Source(a = 0.45,b = 0.04,R = .65)
    geom = Geom(l1=11.6, l2=2.0, l3=3.)
    icg = ICG(source, 5., geom, 0.)
    t = np.arange(-100, 100, 1.)
    from matplotlib import pyplot as plt
    f1,ax1 = plt.subplots()
    ax1.plot(t, icg.ICG(t))
    ax1.set_xlabel('t')
    f1.show()
    return

def test2():
    E = np.arange(-10, 90, .2)
    geom = Geom(l1=11.6, l2=2.0, l3=3.)
    y = resolution(E, 100., 30., 0.45, 0.04, .65, sigma=5., t0=3.5, geom=geom)
    # print E[y!=y]
    from matplotlib import pyplot as plt
    plt.figure()
    plt.plot(E, y)
    plt.show()
    return

def main():
    test1()
    test2()
    return

if __name__ == '__main__': main()

# End of file
