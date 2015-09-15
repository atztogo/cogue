import unittest

from cogue.task import TaskBase, TaskElement, TaskSet

class Test__INIT__(unittest.TestCase):

    def setUp(self):
        pass
    
    def tearDown(self):
        pass

    def test_TaskBase(self):
        self._task = TaskBase()
        print self._task
    
    def test_TaskElement(self):
        self._task = TaskElement()
        print self._task
    
    def test_TaskSet(self):
        self._task = TaskSet(directory="hoge",
                             name="moge")
        print self._task

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(Test__INIT__)
    unittest.TextTestRunner(verbosity=2).run(suite)
    # unittest.main()
