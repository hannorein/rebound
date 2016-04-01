import numpy as np

class Mcmc(object):
    def __init__(self, initial_state, obs):
        self.state = initial_state
        self.obs = obs

    def step(self):
        pass    

class Mh(Mcmc):
    def __init__(self, initial_state, obs, scales=None):
        super(Mh,self).__init__(initial_state, obs)
        if scales is None:
            self.scales = np.ones(self.state.Nvars)
        else:
            self.scales = scales

    def generate_proposal(self):
        prop = self.state.deepcopy()
        shift = self.scales*np.random.normal(size=self.state.Nvars)
        prop.shift_params(shift)
        return prop

    def step(self):
        logp = self.state.get_logp(self.obs)
        proposal = self.generate_proposal() 
        logp_proposal = proposal.get_logp(self.obs)
        if np.exp(logp_proposal-logp)>np.random.uniform():
            self.state = proposal
            return True
        return False

    def step_force(self):
        tries = 1
        while self.step()==False:
            tries += 1
            pass
        return tries



