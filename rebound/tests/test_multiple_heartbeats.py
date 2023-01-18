import rebound
import unittest

class TestMultipleHeartbeats(unittest.TestCase):
    def test_add_heartbeat(self):
        sim = rebound.Simulation()
        sim.integrator = "leapfrog"
        sim.dt = 1.0

        sim.add(m=1)
        sim.move_to_com()

        ts = []
        def heartbeat(sim):
            ts.append(sim.t)

        sim.add_heartbeat(heartbeat)

        n_steps = 40
        sim.integrate(n_steps * sim.dt)

        # + 1 from always running at t=0.0 
        self.assertEqual(n_steps+1, len(ts))
        self.assertEqual(0.0, ts[0])
        self.assertEqual(sim.t, ts[-1])

    def test_add_heartbeat_interval(self):
        sim = rebound.Simulation()
        sim.integrator = "leapfrog"
        sim.dt = 1.0

        sim.add(m=1)
        sim.move_to_com()

        ts = []
        def heartbeat(sim):
            ts.append(sim.t)

        hb_interval = 7.27
        sim.add_heartbeat(heartbeat, interval=hb_interval)

        n_steps = 40
        sim.integrate(n_steps * sim.dt)

        # + 1 from always running at t=0.0 
        self.assertEqual(sim.t // hb_interval + 1, len(ts))
        self.assertEqual(0.0, ts[0])

    def test_add_many_heartbeats_with_phase(self):
        sim = rebound.Simulation()
        sim.integrator = "leapfrog"
        sim.dt = 1.0

        sim.add(m=1)
        sim.move_to_com()

        ts = []
        def heartbeat(sim):
            ts.append(sim.t)

        hb_N = 256
        hb_interval = 256 * sim.dt

        for i in range(hb_N):
            sim.add_heartbeat(heartbeat, interval=hb_interval, phase=i/hb_N)

        n_steps = 300
        sim.integrate(n_steps * sim.dt)

        # + hb_N because all run at 0.0
        self.assertEqual(hb_N + n_steps, len(ts))
        self.assertEqual(0.0, ts[0])
        self.assertEqual(sim.t, ts[-1])


if __name__ == "__main__":
    unittest.main()
