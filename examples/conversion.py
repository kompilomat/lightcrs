from lightcrs import LatLon
from lightcrs import MGRS
from lightcrs import UTM 

if __name__ == "__main__":
    coords = LatLon(48.210033, 16.363449)    # Vienna
    utm_coords = coords.to_UTM()
    mgrs_coords = coords.to_MGRS()

    print(coords)
    print(utm_coords)
    print(mgrs_coords)
    
    mc = MGRS.parse("33UXP012405")
    print(mc.to_LatLon())

    uc = UTM(33, "U", 602066.75288212, 5340350.28679161)
    print(uc.to_LatLon(), uc.to_MGRS())