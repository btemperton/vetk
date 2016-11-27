from unittest import TestCase

import ViralEcologyToolkit


class TestTest(TestCase):
    def testIsString(self):
        s = ViralEcologyToolkit.test()
        self.assertTrue(isinstance(s, str))


class TestCheckFile(TestCase):
    def testNoFile(self):
        with self.assertRaises(IOError):
            ViralEcologyToolkit.checkFasta('foo')
