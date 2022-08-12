"""Test of structure optimization element task."""
import unittest

import cogue.calculator.vasp as vasp


class TestVasp(unittest.TestCase):
    """Test of structure optimization element task."""

    def setUp(self):
        """Set up."""
        self._task = vasp.StructureOptimizationElement()

    def tearDown(self):
        """Tear down."""
        pass

    def test_collect(self):
        """Test."""
        for t in self._task:
            pass
        print(self._task.get_current_cell())


if __name__ == "__main__":
    unittest.main()
