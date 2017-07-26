import unittest

from PysamTestModule_link_pre_489 import build_read

        
class TestModule(unittest.TestCase):

    def test_pass_if_module_can_be_called(self):
        read = build_read()
        self.assertEqual(read.query_name, "hello")
        self.assertEqual(read.query_sequence, "ACGT")
        

if __name__ == "__main__":
    unittest.main()
