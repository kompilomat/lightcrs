import math
import re
from typing import TypeVar, Type

from .constants import *
import lightcrs.latlon as latlon
import lightcrs.utm as utm


# Inline documentation
#
# An MGRS datum is structured as follows
# example from https://en.wikipedia.org/wiki/Military_Grid_Reference_System
#    04Q FJ 12345 67890
#    ---                grid zone designator (GZD) 6x8 degrees
#        --             grid square id, precision level 100 km 
#           -----       easting
#                 ----- northing
#
#    easting/northing precision levels based on digits
#    1 1          10 km square    (precision 1)
#    12 12        1 km square     (precision 2)
#    123 123      100 m square    (precision 3)
#    1234 1234    10 m square     (precision 4)
#    12345 12345  1 m square      (precision 5) 
#
# TRUNCATE don't ROUND
# 


# workaround for typing a classmethod factory
M = TypeVar('M', bound="MGRS")

class MGRS(object):
    def __init__(self,
                gzd : str,
                square_id : str,
                easting : int,
                northing : int,
                precision : int = 5) -> None:

        if len(gzd) < 3:
            self.gzd = "{:0>3}".format(gzd)
        else:
            self.gzd = gzd

        self.square_id = square_id
        self.easting = easting
        self.northing = northing
        self.precision = precision
        self._hash = hash((gzd, square_id, self.easting, self.northing))

    def __hash__(self) -> int:
        return self._hash

    def __eq__(self, other) -> bool:
        if isinstance(other, MGRS):
            return self.gzd == other.gzd and self.square_id == other.square_id \
               and self.easting == other.easting and self.northing == other.northing \
               and self.precision == other.precision
        else:
            return False


    def __repr__(self) -> str:
        if self.precision != 5:
            return f"MGRS('{self.gzd}', '{self.square_id}', {self.easting}, {self.northing}, {self.precision})"
        else:
            return f"MGRS('{self.gzd}', '{self.square_id}', {self.easting}, {self.northing})"

    def __str__(self) -> str:
        easting = self.easting
        northing = self.northing
        for _ in range(self.precision, 5):
            easting = easting // 10
            northing = northing // 10
        easting = str(easting)
        northing = str(northing)
        return "".join((self.gzd, self.square_id, \
                        "0"*(self.precision - len(easting)), easting, 
                        "0"*(self.precision - len(northing)), northing))
   
        
    # returns a tuple in (meters east, meters north)
    def distance(self, other) -> tuple:
        if isinstance(other, MGRS) and self.precision == other.precision:
            factor = 10**(self.precision - 5)
            if self.gzd == other.gzd and self.square_id == other.square_id:
                east = (other.easting - self.easting) * factor
                north = (other.northing - self.northing) * factor
                return east, north
            else:
                raise RuntimeError("Cannot compute distance across squares yet")

        else:
            raise RuntimeError(f"Cannot compare with {other}")

    def to_UTM(self):
        band = self.gzd[-1]
        zone = int(self.gzd[:-1])
        square_e, square_n = self.square_id

        hemisphere = "N" if band >= "N" else "S"

        column = easting_100k_letters[(zone-1) % 3].index(square_e) + 1
        e100k = column * 100e3

        row = northing_100k_Letters[(zone-1) % 2].index(square_n)
        n100k = row * 100e3

        band_latitude = (mgrs_lat_bands.index(band) - 10) * 8
        band_utm = latlon.LatLon(band_latitude, 3).to_UTM()
        band_bottom_northing = math.floor(band_utm.northing / 100e3) * 100e3

        offset = 0
        while (n100k + self.northing + offset < band_bottom_northing):
            offset += 2000e3    # 2,000 km blocks

        return utm.UTM(zone, hemisphere, e100k + self.easting, n100k + offset + self.northing)


    def to_LatLon(self):
        return self.to_UTM().to_LatLon()




    _prefix = re.compile("^(\d?\d[A-Z]) ?([A-Z]{2})")    

    @classmethod
    def parse(cls, strng : str) -> M:
        match = cls._prefix.match(strng)
        if match:
            gzd, square_id = match.group(1, 2)
            offset = len(gzd) + len(square_id)
            remainder = strng[offset:].lstrip().rstrip()
            if  " " in remainder:
                easting, northing = remainder.split(" ")
            else:
                split = len(remainder) // 2
                easting, northing = remainder[:split], remainder[split:]
            
            if len(easting) != len(northing) or len(easting) > 5 or len(northing) > 5:
                raise RuntimeError(f"MGRS easting {easting} and northing {northing} have bad lengths")

            precision = len(easting)
            easting = int(easting)
            northing = int(northing)
            for _ in range (precision, 5):
                easting *= 10
                northing *= 10

            return MGRS(gzd, square_id, easting, northing, precision)
        else:
            raise RuntimeError("Invalid MGRS string")







