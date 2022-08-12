"""Test Task base classes."""
import unittest

from cogue.task import TaskBase, TaskElement, TaskSet


class Test__INIT__(unittest.TestCase):
    """Test TaskBase class."""

    def setUp(self):
        """Set up."""
        pass

    def tearDown(self):
        """Tear down."""
        pass

    def test_TaskBase(self):
        """Test TaskBase class."""
        self._task = TaskBase()
        print(self._task)

    def test_TaskElement(self):
        """Test TaskElement class."""
        self._task = TaskElement()
        print(self._task)

    def test_TaskSet(self):
        """Test TaskSet class."""
        self._task = TaskSet(directory="hoge", name="moge")
        print(self._task)


if __name__ == "__main__":
    suite = unittest.TestLoader().loadTestsFromTestCase(Test__INIT__)
    unittest.TextTestRunner(verbosity=2).run(suite)
    # unittest.main()
