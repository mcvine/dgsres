
from mcvine.workflow import singlextal as sx

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

    
