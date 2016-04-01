import numpy as np

class Mcmc(object):
    def __init__(self, initial_state, obs):
        self.state = initial_state
        self.obs = obs

    def update_logp(self):
        self.logp = -state.get_chi2(self.obs)

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
        prop.shift(shift)
        return prop


    def step(self):
        if self.logp is None:
            self.update_logp()
        proposal = self.generate_proposal() 


