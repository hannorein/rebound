import rebound
import unittest
import ctypes
from rebound.binarydata_field_descriptor import BinarydataFieldDescriptor


class TestBinaryFields(unittest.TestCase):
    def test_unique(self):
        fd_pointer_t = ctypes.POINTER(BinarydataFieldDescriptor)
        fd_pointer = (BinarydataFieldDescriptor*3).in_dll(rebound.clibrebound, "reb_binarydata_field_descriptor_list")
        fd_pointer = ctypes.cast(fd_pointer, fd_pointer_t) # not sure why I have to do it this way
        l = []
        i=0
        while True:
            l.append(fd_pointer[i])
            if fd_pointer[i].name == b'':
                break
            i += 1
        names = []
        for i in l:
            names.append(i.name)
        self.assertEqual(len(names), len(set(names)))
        self.assertGreater(len(names), 50)


if __name__ == "__main__":
    unittest.main()

