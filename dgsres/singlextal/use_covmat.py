# coding: utf-8

# # Resolution Calculation using Covariance Matrix

# from matplotlib import pyplot as plt
import numpy as np, mcvine
import mcvine.cli
# from mcvine_workflow.DGS import ARCS
# import histogram.hdf as hh, histogram as H


# ## Formula
# See NIMA 736(2014)31-39


from . import instrument, pixel
class tofwidths:

    def __init__(self, P, M):
        self.P = P # P chopper
        self.M = M # M chopper
        return

class beamdivs:
    "beam divergence along theta and phi"
    def __init__(self, theta, phi):
        self.theta = theta
        self.phi = phi
        return

def compute(
        sample_yml, Ei, dynamics, psi_scan, instrument, pixel, tofwidths, beamdivs, samplethickness,
        plot=False):
    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    # should P be the T0 chopper?
    L_PM=mcvine.units.parse(instrument.L_m2fc)/mcvine.units.meter # P chopper to M chopper distance
    L_PS= mcvine.units.parse(instrument.L_m2s)/mcvine.units.meter  # P chopper to sample
    L_MS=L_PS-L_PM
    #
    R = mcvine.units.parse(instrument.detsys_radius)/mcvine.units.meter # 

    hkl0 = dynamics.hkl0
    hkl_dir = dynamics.hkl_dir # projection
    psimin = psi_scan.min
    psimax = psi_scan.max
    dpsi = psi_scan.step

    # dynamics calculations
    E = dynamics.E
    dq = dynamics.dq
    hkl = hkl0 + dq*hkl_dir

    from mcni.utils import conversion as Conv
    vi = Conv.e2v(Ei)

    ti = L_PM/vi*1e6 # microsecond
    Ef = Ei - E
    vf = Conv.e2v(Ef)

    # find the psi angle
    from mcvine_workflow.singlextal.io import loadXtalOriFromSampleYml
    xtalori = loadXtalOriFromSampleYml(sample_yml)
    from mcvine_workflow.singlextal.solve_psi import solve
    results = solve(
        xtalori, Ei, hkl, E, psimin, psimax,
        Nsegments = 100)
    from mcvine_workflow.singlextal.coords_transform import hkl2Q
    for r in results:
        xtalori.psi = r*np.pi/180
        print "psi=%s, Q=%s" % (r, hkl2Q(hkl, xtalori))
        print "hkl2Q=%r\n(Q = hkl dot hkl2Q)" % (xtalori.hkl2cartesian_mat(),)
    # these are the psi angles that the particular point of interest will be measured
    # print results
    assert len(results)
    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    # only select the first one. this is OK for most cases but there are cases where more than
    # one psi angles satisfy the condition
    psi = results[0]
    xtalori.psi = psi*np.pi/180
    Q = hkl2Q(hkl, xtalori)
    hkl2Q_mat = xtalori.hkl2cartesian_mat()
    # print Q
    # print hkl2Q_mat
    #
    Q_len = np.linalg.norm(Q); print Q_len
    ki = Conv.e2k(Ei); print ki
    kiv = np.array([ki, 0, 0])
    kfv = kiv - Q; print kfv
    #
    # ** Verify the momentum and energy transfers **
    # print Ei-Conv.k2e(np.linalg.norm(kfv))
    # print Ei-Ef
    assert np.isclose(Ei-Ef, E)

    # ** Compute detector pixel position **
    z = kfv[2]/(kfv[0]**2+kfv[1]**2)**.5 * R
    L_SD=(z**2+R**2)**.5
    # print z, L_SD

    # ### Constants
    eV = 1.60218e-19
    meV = eV*1e-3
    mus = 1.e-6
    hbar= 1.0545718e-34
    AA = 1e-10
    m = 1.6750e-24 * 1e-3 #kg
    from numpy import sin, cos

    # dE calcuation starts here
    # ## Differentials
    pE_pt = -m*(vi**3/L_PM + vf**3/L_SD * L_MS/L_PM)
    # convert to eV/microsecond
    pE_pt /= meV/mus
    # print pE_pt

    pE_ptMD = m*vf**3/L_SD
    pE_ptMD /= meV/mus
    # print pE_ptMD

    pE_pLPM = m/L_PM * (vi**2 + vf**3/vi * L_MS/ L_SD)
    pE_pLPM /= meV
    # print pE_pLPM

    pE_pLMS= -m/L_SD * (vf**3/vi)
    pE_pLMS /= meV
    # print pE_pLMS

    pE_pLSD = -m*vf*vf/L_SD
    pE_pLSD /= meV
    # print pE_pLSD

    # we don't need pE_pLSD, instead we need pE_pR and pE_pz. R and z are cylinder radius and z coordinate
    pE_pR = pE_pLSD * (R/L_SD)
    pE_pz = pE_pLSD * (z/L_SD)
    # print pE_pR, pE_pz

    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    # ## ** Paramters: Estimate of standard deviations
    # tau_P = 10 # microsecond
    # tau_M = 8 # microsecond
    tau_P = tofwidths.P
    tau_M = tofwidths.M
    #
    # tau_D = 10 # microsecond
    tau_D = mcvine.units.parse(pixel.radius)/mcvine.units.meter*2/vf*1e6 # microsecond

    # ## Calculations
    pE_p_vec = [pE_pt, pE_ptMD, pE_pLPM, pE_pLMS, pE_pR, pE_pz]
    pE_p_vec = np.array(pE_p_vec)
    J_E = pE_p_vec/E

    # print J_E
    sigma_t = (tau_P**2+tau_M**2)**.5
    sigma_tMD = (tau_M**2+tau_D**2)**.5
    # XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    div = (beamdivs.theta**2 + beamdivs.phi**2)**.5  # a crude approx
    sigma_LPM = L_PM * div*div

    # mainly due to sample size
    sigma_LMS = samplethickness

    # mainly due to det tube diameter
    sigma_R = mcvine.units.parse(pixel.radius)/mcvine.units.meter*2

    # pixel size
    sigma_z = mcvine.units.parse(pixel.height)/mcvine.units.meter

    sigma = np.array([sigma_t, sigma_tMD, sigma_LPM, sigma_LMS, sigma_R, sigma_z])
    # print sigma
    sigma2 = sigma*sigma
    # print sigma2
    sigma2 = np.diag(sigma2)

    # print J_E
    # print np.dot(sigma2, J_E)
    cov = np.dot(J_E, np.dot(sigma2, J_E))

    # print cov, np.sqrt(cov)
    sigma_E = E*np.sqrt(cov)
    # print sigma_E

    # Not sure if this is right:
    # 
    # * ** Note: this may be more like FWHM than sigma_E because of the approx I made **
    # * ** FWHM is 2.355 sigma **

    # ## Include Q
    print "ti=",ti
    tf = L_SD/vf*1e6; print "tf=",tf
    thetai = 0
    phii = 0
    print "R=", R
    print "Q=", Q

    eeta = np.arctan2(Q[1], Q[0])
    print "eeta=", eeta

    pQx_pt = -m/hbar*(L_PM/ti/ti/mus/mus*cos(thetai)*cos(phii)
                      +R/tf/tf/mus/mus*L_MS/L_PM*cos(eeta))
    pQx_pt/=1./AA/mus
    # print pQx_pt

    pQx_ptMD = m/hbar * R/tf/tf * cos(eeta) / mus/mus
    pQx_ptMD /= 1./AA/mus
    # print pQx_ptMD

    pQx_pLPM = m/hbar *(cos(thetai) * cos(phii)/ti + ti/tf/tf * R*L_MS/L_PM/L_PM * cos(eeta)) / mus
    pQx_pLPM /= 1./AA
    # print pQx_pLPM

    pQx_pLMS = -m/hbar * R/tf/tf*ti/L_PM*cos(eeta) / mus
    pQx_pLMS /= 1./AA
    # print pQx_pLMS

    pQx_pR = -m/hbar/tf*cos(eeta) / mus
    pQx_pR /= 1./AA
    # print pQx_pR

    pQx_peeta = m/hbar * R/tf*sin(eeta) /mus
    pQx_peeta /= 1./AA
    # print pQx_peeta

    pQx_pthetai = -m/hbar*L_PM/ti*sin(thetai)*cos(phii)/mus
    pQx_pthetai/=1./AA
    # print pQx_pthetai

    pQx_pphii = -m/hbar*L_PM/ti*cos(thetai)*sin(phii)/mus
    pQx_pphii/=1./AA
    # print pQx_pphii

    pQx_p_vec = [pQx_pt, pQx_ptMD, pQx_pLPM, pQx_pLMS, pQx_pR, 0, pQx_peeta, pQx_pthetai, pQx_pphii]
    pQx_p_vec = np.array(pQx_p_vec)
    J_Qx = pQx_p_vec/Q_len

    # **Qy**
    pQy_pt = -m/hbar*(L_PM/ti/ti*sin(thetai)*cos(phii)+R/tf/tf*L_MS/L_PM*sin(eeta))/mus/mus
    pQy_pt/=1./AA/mus
    # print pQy_pt

    pQy_ptMD = m/hbar * R/tf/tf * sin(eeta) / mus/mus
    pQy_ptMD /= 1./AA/mus
    # print pQy_ptMD

    pQy_pLPM = m/hbar *(sin(thetai) * cos(phii)/ti 
                        + ti/tf/tf * R*L_MS/L_PM/L_PM * sin(eeta)) / mus
    pQy_pLPM /= 1./AA
    # print pQy_pLPM

    pQy_pLMS = -m/hbar * R/tf/tf*ti/L_PM*sin(eeta) / mus
    pQy_pLMS /= 1./AA
    # print pQy_pLMS

    pQy_pR = -m/hbar/tf*sin(eeta) / mus
    pQy_pR /= 1./AA
    # print pQy_pR

    pQy_peeta = -m/hbar * R/tf*cos(eeta) /mus
    pQy_peeta /= 1./AA
    # print pQy_peeta

    pQy_pthetai = m/hbar*L_PM/ti*cos(thetai)*cos(phii)/mus
    pQy_pthetai/=1./AA
    # print pQy_pthetai

    pQy_pphii = -m/hbar*L_PM/ti*sin(thetai)*sin(phii)/mus
    pQy_pphii/=1./AA
    # print pQy_pphii

    pQy_p_vec = [pQy_pt, pQy_ptMD, pQy_pLPM, pQy_pLMS, pQy_pR, 0, pQy_peeta, pQy_pthetai, pQy_pphii]
    pQy_p_vec = np.array(pQy_p_vec)
    J_Qy = pQy_p_vec/Q_len

    # ** Qz **
    pQz_pt = -m/hbar*(L_PM/ti/ti*sin(phii)+z/tf/tf*L_MS/L_PM)/mus/mus
    pQz_pt/=1./AA/mus
    # print pQz_pt

    pQz_ptMD = m/hbar * z/tf/tf /mus/mus
    pQz_ptMD /= 1./AA/mus
    # print pQz_ptMD

    pQz_pLPM = m/hbar *(sin(phii)/ti + ti/tf/tf * z*L_MS/L_PM/L_PM) / mus
    pQz_pLPM /= 1./AA
    # print pQz_pLPM

    pQz_pLMS = -m/hbar * z/tf/tf*ti/L_PM / mus
    pQz_pLMS /= 1./AA
    # print pQz_pLMS

    pQz_pz = -m/hbar/tf / mus
    pQz_pz/=1./AA
    # print pQz_pz

    pQz_pphii = m/hbar*L_PM/ti*cos(phii)/mus
    pQz_pphii/=1./AA
    # print pQz_pphii

    pQz_p_vec = [pQz_pt, pQz_ptMD, pQz_pLPM, pQz_pLMS, 0, pQz_pz, 0, 0, pQz_pphii]
    pQz_p_vec = np.array(pQz_p_vec)
    J_Qz = pQz_p_vec/Q_len

    # ** Here we need to extend the J vector for E to include the additional variables eeta, thetai, and phii **
    pE_p_vec = [pE_pt, pE_ptMD, pE_pLPM, pE_pLMS, pE_pR, pE_pz,0,0,0]
    pE_p_vec = np.array(pE_p_vec)
    J_E = pE_p_vec/E
    J = np.array( (J_Qx, J_Qy, J_Qz, J_E) )

    # ## ** Parameters
    sigma_eeta = mcvine.units.parse(pixel.radius)/mcvine.units.parse(instrument.detsys_radius)
    # sigma_thetai = 0.01
    sigma_thetai = beamdivs.theta
    # sigma_phii = 0.01
    sigma_phii = beamdivs.phi
    sigma = np.array([sigma_t, sigma_tMD, sigma_LPM, sigma_LMS, sigma_R, sigma_z, sigma_eeta, sigma_thetai, sigma_phii])
    sigma2 = sigma**2
    sigma2 = np.diag(sigma2)

    # print J.shape, sigma2.shape
    cov = np.dot(J, np.dot(sigma2, J.T))
    # print cov

    M = np.linalg.inv(cov)
    # print M

    np.dot(cov, M)

    # ## Ellipsoid
    # hkl = hkl0+hkl_dir*x
    # dh,dk,dl = dx2dhkl*dx 
    dx2dhkl = np.array(hkl_dir)

    # dQ = dx * dx2dhkl dot hkl2Q
    # so dx2dQ = dx2dhkl * hkl2Q
    # dQ = dx * dx2dQ
    dx2dQ = np.dot(dx2dhkl, hkl2Q_mat)
    # print dx2dQ

    # [dQx,dQy,dQz,dE] = [dx dE] dot dxdE2dQdE
    L=dxdE2dQdE = np.array([list(dx2dQ) + [0],  
                            [ 0.        ,  0.,          0.,  1]
                            ])

    np.dot([1,1], dxdE2dQdE)

    # $ [dX1,\; dX2,\; dX3,\; dX4]\; M\; [dX1,\; dX2,\; dX3,\; dX4 ]^T = 2ln(2)$
    # 
    # $ dX_i = \frac{dQ_i}{|Q|}$ for i = 1,2,3
    # 
    # $ dX_4 = \frac{dE}{E}$
    # 
    # Let 
    # $ U = diag\big( \frac{1}{|Q|},\; \frac{1}{|Q|},\; \frac{1}{|Q|},\; 1/E \big) $
    # 
    # $ [dx,\; dE]\; L U MU^TL^T [dx,\; dE ]^T = 2ln(2)$
    # 
    # Let $N=L U MU^TL^T $
    # print Q_len, E
    U = np.diag([1./Q_len, 1./Q_len, 1./Q_len, 1./E])
    N = LUMUTLT = np.dot(L, np.dot(U, np.dot(M, np.dot(U.T, L.T))))
    # print N
    # print 2*np.log(2)
    r = np.linalg.eig(N)
    mR = r[1]; lambdas = r[0]
    # print np.dot(mR, mR.T)
    np.dot(np.dot(mR.T, N), mR)

    # $ u = [dx,\;dE]$
    # 
    # $ u N u^T = 2ln(2)$        .... (1)
    # 
    # Find eigen values ($\lambda_1$, $\lambda_2$) and eigne vectors ($e_1$, $e_2$, column vectors) of N,
    # and let 
    # 
    # $ R = [e_1,\;e_2] $
    # 
    # Then
    # 
    # $ N' = R^T N R = diag([\lambda_1, \lambda_2]) $
    # 
    # or
    # 
    # $ N = R N' R^T $
    # 
    # With $N'$ we can rewrite (1) as
    # 
    # $ u'N'{u'}^T = 2ln2 = \lambda_1 {u'}_1^2 + \lambda_2 {u'}_2^2 $
    # 
    # where
    # 
    # $ u' = u . R $

    # ${u'}_1 = \sqrt{2ln2/\lambda_1}*cos(\theta)$
    # 
    # ${u'}_2 = \sqrt{2ln2/\lambda_2}*sin(\theta)$

    # In[ ]:

    RR = 2*np.log(2)
    theta = np.arange(0, 360, 1.)*np.pi/180
    u1p = np.sqrt(RR/lambdas[0])*np.cos(theta)
    u2p = np.sqrt(RR/lambdas[1])*np.sin(theta)
    up = np.array([u1p, u2p]).T

    # print up.shape
    u = np.dot(up, mR.T)

    if plot:
        from matplotlib import pyplot as plt
        plt.plot(u[:,0], u[:,1], '.')
        # plt.xlim(-.35, .1)
        # plt.ylim(-5., 5.)
        plt.show()
    # ellipsoid coordinates, eigen vectors and eigen values of the scaled inverse covariance
    return u, mR, lambdas


