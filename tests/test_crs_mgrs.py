import unittest
from random import Random, random

from lightcrs import MGRS


class TestMGRS(unittest.TestCase):
    def setUp(self):
        global random
        random = Random(42)

    def test_mgrs_distance(self):
       a = MGRS.parse("04UGJ2345300456")
       b = MGRS.parse("04UGJ2345500459")
       self.assertEqual(a.distance(b), (2, 3))


if __name__ == '__main__':
    unittest.main()

