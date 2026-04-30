import rebound
import rebound.binary_field_descriptor
import unittest

class TestBinaryFields(unittest.TestCase):
    def test_unique(self):
        l = rebound.binary_field_descriptor.binary_field_descriptor_list()
        names = []
        for i in l:
            names.append(i.name)
        self.assertEqual(len(names), len(set(names)))


if __name__ == "__main__":
    unittest.main()

