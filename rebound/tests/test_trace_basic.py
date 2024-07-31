import rebound
import unittest


class TestIntegratorTraceBasic(unittest.TestCase):
   
    def test_trace_encounter_condition(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=1e-3,a=1)
        sim.add(m=1e-3,a=1.1,f=0.25) # First step will be rejected.
        sim.add(m=1e-3,a=3.1)
        sim.integrator = "TRACE"
        sim.ri_trace.r_crit_hill *= 1.21
        sim.dt = 0.14 

        for k in range(3):
            sim.step()
            self.assertEqual(sim.ri_trace._current_C,0)
            self.assertEqual(sim.ri_trace._encounter_N,3)
            self.assertEqual(sim.ri_trace._encounter_N_active,3)
            for i in range(sim.N):
                for j in range(i+1,sim.N):
                    if i==1 and j==2: 
                        self.assertEqual(sim.ri_trace._current_Ks[i*sim.N+j],1)
                    else:
                        self.assertEqual(sim.ri_trace._current_Ks[i*sim.N+j],0)
        
        # Change condition. No more encounters.
        sim.ri_trace.r_crit_hill = 0.1
        for k in range(3):
            sim.step()
            self.assertEqual(sim.ri_trace._current_C,0)
            self.assertEqual(sim.ri_trace._encounter_N,1)
            for i in range(sim.N):
                for j in range(i+1,sim.N):
                    self.assertEqual(sim.ri_trace._current_Ks[i*sim.N+j],0)

    def test_trace_encounter_prediction(self):
        sim = rebound.Simulation()
        sim.add(m=1)
        sim.add(m=9.55e-4,x=5.2)
        sim.add(x=5.3,y=0.36,vy=-7.2) # Non-encounter prediction misses this
        sim.integrator = "TRACE"
        sim.dt = 0.01
        sim.ri_trace.r_crit_hill = 1
        sim.step()
        
        self.assertEqual(sim.ri_trace._encounter_N,3)



if __name__ == "__main__":
    unittest.main()
