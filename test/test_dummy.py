"""Just a dummy test function that always passes to make sure that the test suite is working."""

import unittest

class TestDummy(unittest.TestCase):
    def test_dummy(self):
        self.assertEqual(1, 1)