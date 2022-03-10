import re
import math

from .constants import *
import lightcrs.utm as utm
import lightcrs.mgrs as mgrs

# WGS84 geoid is assumed
# N .. 0
# E .. 90
# S .. 180
# W .. 270
# latitude south..north [-90..90]
# longitude west..east [-180..180]

class LatLon(object):
    def __init__(self,
                lat : float,
                lon : float) -> None:

        self.lat = lat
        self.lon = lon
        self._hash = hash((self.lat, self.lon))


    def __hash__(self) -> int:
        return self._hash

    def __repr__(self) -> str:
        return f"LatLon({self.lat}, {self.lon})"

    def __str__(self) -> str:
        if self.lat < 0: 
            ns = "S"
            lat = abs(self.lat)
        else:
            ns = "N"
            lat = self.lat

        if self.lon < 0: 
            ew = "W"
            lon = abs(self.lon)
        else:
            ew = "E"
            lon = self.lon       
        return f"{lat}{ns}, {lon}{ew}"


    def to_UTM(self) -> utm.UTM:
        """ Computes Universal Transverse Mercator coordinates 
        from WGS84 based latitude longitude coordinates (e.g. GPS).

        Grid zones are 8 degrees latitude. 
        N0 degrees is offset 10 into latitude bands.

        Args:
            latitude (float): [-90 .. 90] degrees
            longitude (float): [-180 .. 180] degrees

        Returns:
            UTM : namedtuple
        """
        if not -80.0 <= self.lat <= 84.0:
            raise RuntimeError(f"latitude {self.lat} outside UTM limits") 
        if self.lon == 180:
            self.lon = -180.
        zone = math.floor((self.lon + 180) / 6) + 1
        lon_central_meridian = math.radians((zone - 1)*6 - 180 + 3)

        lat_band_idx = int(math.floor((self.lat/8) + 10))
        lat_band = mgrs_lat_bands[lat_band_idx]

        # special case Norway
        if zone == 31 and lat_band == "V" and self.lon >= 3.0:
            zone += 1
            lon_central_meridian += math.radians(6)

        # special case Svalbard
        if zone == 32 and lat_band == "X" and self.lon < 9.0:
            zone -= 1
            lon_central_meridian -= math.radians(6)  
        if zone == 32 and lat_band == "X" and self.lon >= 9.0:
            zone += 1
            lon_central_meridian += math.radians(6)  
        if zone == 34 and lat_band == "X" and self.lon < 21.0:
            zone -= 1
            lon_central_meridian -= math.radians(6)  
        if zone == 34 and lat_band == "X" and self.lon >= 21.0:
            zone += 1
            lon_central_meridian += math.radians(6)  
        if zone == 36 and lat_band == "X" and self.lon < 33.0:
            zone -= 1
            lon_central_meridian -= math.radians(6)  
        if zone == 36 and lat_band == "X" and self.lon >= 33.0:
            zone += 1
            lon_central_meridian += math.radians(6)  

        phi = math.radians(self.lat)
        lam = math.radians(self.lon) - lon_central_meridian

        cos_lam = math.cos(lam)
        sin_lam = math.sin(lam)
        tan_lam = math.tan(lam)


        tau = math.tan(phi)
        sigma = math.sinh(eccentricity * math.atanh(eccentricity * tau / math.sqrt(1 + tau**2)))
        tau_prime = tau * math.sqrt(1 + sigma**2) - sigma * math.sqrt(1 + tau**2)

        xi_prime = math.atan2(tau_prime, cos_lam)
        eta_prime = math.asinh(sin_lam / math.sqrt(tau_prime**2 + cos_lam**2))

        xi = xi_prime
        eta  = eta_prime

        for j in range(1, 7):
            xi  += alpha[j] * math.sin(2* j * xi_prime) * math.cosh(2 * j * eta_prime)
            eta += alpha[j] * math.cos(2* j * xi_prime) * math.sinh(2 * j * eta_prime)

        x = scale * A * eta
        y = scale * A * xi

        #  convergence: Karney 2011 Eq 23, 24
        p_prime = 1
        q_prime = 0
        for j in range(1, 7):
            p_prime += 2 * j * alpha[j] * math.cos(2 * j * xi_prime) * math.cosh(2 * j * eta_prime)
            q_prime += 2 * j * alpha[j] * math.sin(2 * j * xi_prime) * math.sinh(2 * j * eta_prime)

        gamma_prime = math.atan(tau_prime / math.sqrt(1 + tau_prime**2) * tan_lam)
        gamma_pprime = math.atan2(q_prime, p_prime)
        gamma = gamma_prime + gamma_pprime

        # scale: Karney 2011 Eq 25
        sin_phi = math.sin(phi)
        k_prime = math.sqrt(1 - eccentricity**2 * sin_phi**2) * math.sqrt(1 + tau**2) / \
            math.sqrt(tau_prime**2 + cos_lam**2)
        k_pprime = (A / semimajor_axis) * math.sqrt(p_prime**2 + q_prime**2)
        k = scale * k_prime * k_pprime

        # shift origin
        x += 500000.0  # false easting
        if y < 0:
            y += 10000000.0  # false northing

        convergence = math.degrees(gamma)
        hemisphere = "N" if self.lat >= 0.0 else "S"

        # "zone", "band", "hemisphere", "easting", "northing"
        return utm.UTM(zone, hemisphere, x, y)




    def to_MGRS(self, precision=5) -> mgrs.MGRS:
        ucoords = self.to_UTM()

        lat_band_idx = int(math.floor((self.lat/8) + 10))
        band = mgrs_lat_bands[lat_band_idx]

        column = math.floor(ucoords.easting / 100e3)
        square_e = easting_100k_letters[(ucoords.zone-1) % 3][column - 1]
        
        row = math.floor(ucoords.northing / 100e3) % 20
        square_n = northing_100k_Letters[(ucoords.zone-1) % 2][row]
    
        easting = int(ucoords.easting % 100e3)
        northing = int(ucoords.northing % 100e3)

        gzd = "{:0>2}{}".format(ucoords.zone, band)
        square_id = square_e + square_n

        return mgrs.MGRS(gzd, square_id, easting, northing, precision)