
from mcvine.workflow import singlextal as sx

class arcs:
    instrument = sx.instrument(
        name = 'ARCS',
        detsys_radius = "3.*meter",
        L_m2s = "13.6*meter",
        offset_sample2beam = "-0.15*meter" # offset from sample to saved beam. don't change this unless you are sure what you are doing
    )

    pixel = sx.pixel(
        radius = "0.5*inch",
        height = "1.*meter/128",
        pressure = "10*atm",
    )

    @classmethod
    def scattering_angle_constraints(cls, theta, phi):
        return ((theta<135.) * (theta>-28)) * (phi<26) * (phi>-27)

class sequoia:
    instrument = sx.instrument(
        name = 'SEQ',
        detsys_radius = "5.5*meter",
        L_m2s = "20.05*meter",
        offset_sample2beam = "-0.15*meter" # offset from sample to saved beam. don't change this unless you are sure what you are doing
    )

    pixel = sx.pixel(
        radius = "0.5*inch",
        height = "1.2*meter/128",
        pressure = "10*atm",
    )

    @classmethod
    def scattering_angle_constraints(cls, theta, phi):
        return ((theta<60.) * (theta>-30)) * (phi<18) * (phi>-18)

class cncs:
    instrument = sx.instrument(
        name = 'CNCS',
        detsys_radius = "3.5*meter",
        L_m2s = "36.2*meter",
        offset_sample2beam = "-0.15*meter" # offset from sample to saved beam. don't change this unless you are sure what you are doing
    )

    pixel = sx.pixel(
        radius = "0.5*inch",
        height = "2.*meter/128",
        pressure = "6.*atm",
    )

    @classmethod
    def scattering_angle_constraints(cls, theta, phi):
        return ((theta<140.) * (theta>-50)) * (phi<16) * (phi>-16)


def chess_pixel_orientation_func(theta, phi):
    import numpy as np
    from mcni.neutron_coordinates_transformers.mcstasRotations import toMatrix
    if np.rad2deg(np.abs(theta))>34:
        m = np.dot(toMatrix(0, np.rad2deg(phi), 0), toMatrix(0, 0, 90.))
    else:
        m = np.dot(toMatrix(np.rad2deg(theta)-90., 0, 0), toMatrix(0, np.rad2deg(phi), 0))
    return m

class chess:

    instrument = sx.instrument(
        name = 'CHESS',
        detsys_radius = "2.5*meter", detsys_shape = 'sphere',
        L_m2s = "31.5*meter",
        offset_sample2beam = "-0.15*meter", # offset from sample to saved beam. don't change this unless you are sure what you are doing
        pixel_orientation_func = chess_pixel_orientation_func
    )

    pixel = sx.pixel(
        radius = "0.5*inch",
        height = "1.5*meter/128",
        pressure = "6.*atm",
    )

    @classmethod
    def scattering_angle_constraints(cls, theta, phi):
        return ((theta<180.) * (theta>-180)) * (phi<90) * (phi>-90)


class amateras:
    instrument = sx.instrument(
        name = 'AMATERAS',
        detsys_radius = "4.0*meter",
        L_m2s = "30.0*meter",
        offset_sample2beam = "-0.15*meter" # offset from sample to saved beam. don't change this unless you are sure what you are doing
    )

    pixel = sx.pixel(
        radius = "0.5*inch",
        height = "0.02425*meter",
        pressure = "10.0*atm",
    )

    @classmethod
    def scattering_angle_constraints(cls, theta, phi):
        return ((theta<140.) * (theta>-40)) * (phi<20.5) * (phi>-20.5)

class LET:
    instrument = sx.instrument(
        name = 'LET',
        detsys_radius = "3.5*meter",
        L_m2s = "25*meter",
        offset_sample2beam = "-0.15*meter" # offset from sample to saved beam. don't change this unless you are sure what you are doing
    )

    pixel = sx.pixel(
        radius = "0.5*inch",
        height = "4./1024*meter",
        pressure = "10.0*atm",
    )

    @classmethod
    def scattering_angle_constraints(cls, theta, phi):
        return ((theta<140.) * (theta>-40.)) * (phi<30.) * (phi>-30.)
