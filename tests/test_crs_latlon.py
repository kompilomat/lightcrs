import unittest
from random import Random, random

from lightcrs import LatLon

# needs to be installed: pip install mgrs
import mgrs



class TestLatLon(unittest.TestCase):
    def setUp(self):
        global random
        random = Random(42)
        self.mgrs_reference = mgrs.MGRS()


    def test_random_mgrs(self):
       for _ in range(2000):
            lat = random.randint(-80, 84)
            lon = random.randint(-180,180)
            goal = self.mgrs_reference.toMGRS(lat, lon)
            impl = LatLon(lat, lon).to_MGRS()
            self.assertEqual(goal, str(impl), f"lat: {lat} lon: {lon}")


if __name__ == '__main__':
    unittest.main()
