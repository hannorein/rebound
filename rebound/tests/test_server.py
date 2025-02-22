import rebound
import unittest
import os
import warnings
import urllib.request


class TestServer(unittest.TestCase):
    def test_start_without_download(self):
        with open("rebound.html", "w") as f:
            f.write("<h1>Hi</h1>")
        sim = rebound.Simulation()
        sim.start_server(port=1234)
        self.assertNotEqual(sim._server_data,None)
        self.assertEqual(sim._server_data.contents.ready,1)
        self.assertEqual(sim._server_data.contents.port,1234)
    
    def test_start_with_download(self):
        if os.path.isfile("rebound.html"):
            os.remove("rebound.html")
        sim = rebound.Simulation()
        with warnings.catch_warnings(record=True) as w: 
            warnings.simplefilter("always")
            sim.start_server(port=1234)
            self.assertGreater(len(w),0)
        self.assertNotEqual(sim._server_data,None)
        self.assertEqual(sim._server_data.contents.ready,1)
        self.assertEqual(sim._server_data.contents.port,1234)
    
    def test_connect(self):
        teststring = "<h1>Hi</h1>"
        with open("rebound.html", "w") as f:
            f.write(teststring )
        sim = rebound.Simulation()
        sim.start_server(port=1234)
        contents = urllib.request.urlopen("http://localhost:1234/").read()
        contents = contents.decode("ascii")
        self.assertEqual(contents, teststring)

    def test_pause(self):
        teststring = "<h1>Hi</h1>"
        with open("rebound.html", "w") as f:
            f.write(teststring )
        sim = rebound.Simulation()
        sim.start_server(port=1234)
        sim._status = -1;
        contents = urllib.request.urlopen("http://localhost:1234/keyboard/32")
        self.assertEqual(sim._status, -3)

    
    
if __name__ == "__main__":
    unittest.main()
