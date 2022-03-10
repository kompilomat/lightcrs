import math

from .constants import *

import lightcrs.latlon as latlon
import lightcrs.mgrs as mgrs

class UTM(object):
    def __init__(self,
                zone : int,
                hemisphere : str,
                easting : float,
                northing : float) -> None:
        
        self.zone = zone
        self.hemisphere = hemisphere
        self.easting = easting
        self.northing = northing
        self._hash = hash((self.zone, self.hemisphere, self.easting, self.northing))


    def __hash__(self) -> int:
        return self._hash

    def __repr__(self) -> str:
        return f"UTM({self.zone}, '{self.hemisphere}', {self.easting}, {self.northing})"

    def __str__(self) -> str:
        return f"{self.zone} {self.hemisphere} {self.easting} {self.northing}"


    def epsg(self) -> int:
        if self.hemisphere == "N":
            return 32600 + self.zone
        elif self.hemisphere == "S":
            return 32700  + self.zone


    def to_LatLon(self):
        false_easting = 500e3
        false_northing = 10000e3

        # x, y relative central meridian and equator
        x = self.easting - false_easting
        if self.hemisphere == "S":
            y = self.northing - false_northing
        else:
            y = self.northing

        #  from Karney 2011 Eq 15-22, 36
        eta = x / (scale * A)
        xi = y / (scale * A)

        xi_prime = xi
        eta_prime = eta
        for j in range(1, 7):
            xi_prime -= beta[j] * math.sin(2*j*xi) * math.cosh(2*j*eta)
            eta_prime -= beta[j] * math.cos(2*j*xi) * math.sinh(2*j*eta)

        sinh_eta_prime = math.sinh(eta_prime)
        sin_xi_prime = math.sin(xi_prime)
        cos_xi_prime = math.cos(xi_prime)

        tau_prime = sin_xi_prime / math.sqrt(sinh_eta_prime**2 + cos_xi_prime**2)


        tau = tau_prime

        # first differential update
        sigma_i = math.sinh(eccentricity * math.atanh(eccentricity*tau / math.sqrt(1 + tau**2)))
        tau_i_prime = tau * math.sqrt(1 + sigma_i**2) - sigma_i * math.sqrt(1 + tau**2)
        dtau = (tau_prime - tau_i_prime) / math.sqrt(1 + tau_i_prime**2) \
             * (1 + (1 - eccentricity**2) * tau**2) \
             / ((1 - eccentricity**2) * math.sqrt(1 + tau**2))
        tau += dtau

        while abs(dtau) > 1e-12:
            sigma_i = math.sinh(eccentricity * math.atanh(eccentricity*tau / math.sqrt(1 + tau**2)))
            tau_i_prime = tau * math.sqrt(1 + sigma_i**2) - sigma_i * math.sqrt(1 + tau**2)
            dtau = (tau_prime - tau_i_prime) / math.sqrt(1 + tau_i_prime**2) \
                 * (1 + (1 - eccentricity**2) * tau**2) \
                 / ((1 - eccentricity**2) * math.sqrt(1 + tau**2))
            tau += dtau

        # tau should have converged now
        phi = math.atan(tau)
        lbda = math.atan2(sinh_eta_prime, cos_xi_prime)

        # convergence: Karney 2011 Eq 26, 27
        p = 1
        q = 0
        for j in range(1, 7):
            p -= 2*j*beta[j] * math.cos(2*j*xi) * math.cosh(2*j*eta)
            q += 2*j*beta[j] * math.sin(2*j*xi) * math.sinh(2*j*eta)

        gamma_prime = math.atan(math.tan(xi_prime) * math.tanh(eta_prime))
        gamma_pprime = math.atan2(q, p)
        gamma = gamma_prime + gamma_pprime

        # scale: Karney 2011 Eq 28
        sin_phi = math.sin(phi)
        k_prime = math.sqrt(1 - eccentricity**2 * sin_phi**2) \
                * math.sqrt(1 + tau**2) \
                * math.sqrt(sinh_eta_prime**2 + cos_xi_prime**2)
        k_pprime = A / semimajor_axis / math.sqrt(p**2 + q**2)
        k = scale * k_prime * k_pprime
        # k is scale

        lbda += math.radians((self.zone - 1)*6 - 180 + 3)

        lat = math.degrees(phi)
        lon = math.degrees(lbda)
        convergence = math.degrees(gamma)

        return latlon.LatLon(lat, lon)
        


    def to_MGRS(self, precision=5) -> mgrs.MGRS:
        coords = self.to_LatLon()

        lat_band_idx = int(math.floor((coords.lat/8) + 10))
        band = mgrs_lat_bands[lat_band_idx]

        column = math.floor(self.easting / 100e3)
        square_e = easting_100k_letters[(self.zone-1) % 3][column - 1]
        
        row = math.floor(self.northing / 100e3) % 20
        square_n = northing_100k_Letters[(self.zone-1) % 2][row]
    
        easting = int(self.easting % 100e3)
        northing = int(self.northing % 100e3)

        gzd = "{:0>2}{}".format(self.zone, band)
        square_id = square_e + square_n

        return mgrs.MGRS(gzd, square_id, easting, northing, precision)
