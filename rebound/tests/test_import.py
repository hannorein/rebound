import warnings
import unittest

class TestImport(unittest.TestCase):
    def test_no_warning_on_import(self):
        with warnings.catch_warnings(record=True) as w: 
            warnings.simplefilter("always")
            import rebound
            self.assertEqual(0,len(w))
    
if __name__ == "__main__":
    unittest.main()
