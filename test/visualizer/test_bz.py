import unittest

import numpy as np
from mayavi import mlab

from cogue.visualizer.bz import VisualizeBrillouinZone


class TestBz(unittest.TestCase):
    def setUp(self):
        lattice = [[4.65, 0, 0], [0, 4.75, 0], [0, 0, 3.25]]
        self._bz = VisualizeBrillouinZone(np.linalg.inv(lattice))

    def tearDown(self):
        pass

    def test_bz(self):
        self.assertTrue(self._bz.run(with_axis=True))


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(TestBz)
    unittest.TextTestRunner(verbosity=2).run(suite)
    # unittest.main()
    mlab.show()
