P0 = 101325.0
T0 = 288.15
G0 = 9.80665
R  = 287.04
GAMMA = 1.4

class FlightConditions(object):
    def __init__(self, speed, altitude, dT=0.0):
        self.g = G0
        self.atm = ISAtmosphere(altitude, dT)
        self.velocity = float(speed)
        self.Mach = self.velocity / self.atm.sound_speed
        self.dynamic_pressure = 0.5 * self.atm.density * self.velocity**2.0
        
class ISAtmosphere:
    def __init__(self, altitude, dT=0):
        altitude = float(altitude)
        T = T0 - 6.5 * altitude / 1000 + dT
        P = P0 * (1 - 0.0065 * altitude / T0) ** 5.2561
        rho = P / (R * T)
        a = (GAMMA * R * T) ** 0.5
        self.dT          = dT
        self.altitude    = altitude
        self.temperature = T
        self.pressure    = P
        self.sound_speed = a
        self.density     = rho