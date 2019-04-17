
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
        L_m2fc = "18*meter",
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

    class violini:

        tau_P = 10
        tau_M = 8
        sigma_thetai = 0.01
        sigma_phii = 0.01

    
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

    
