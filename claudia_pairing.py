'''
Basic model neuron with a synapse subjected to STDP.
Spike dynamics follow HH.
STDP rule outlined in Clopath et al., Nat. Neurosci. 13(3), 344-352 (2010) 
Protocol reproducing plasticity data from Sjoestroem et al. Neuron 32, 1149-1164 (2001); 15 times 5 pairings repeated at 0.1, 10, 20, 40, 50 Hz, delata t = +10ms, -10ms

Ben Torben-Nielsen, Hebrew University
Claudia Clopath, somewhere in NY
'''

import numpy as np
import matplotlib.pyplot as plt
import neuron
from neuron import h

NO_REPS = 5 # 5 pairings. fixed because the result is multiplied by 15 (15x5=75 as in the experiments)
DT=0.1 # ms, set the integration step. Important for the DELAY_STEPS in stdp_cc.mod
POST_AMP = 10 # nA, amplitude of current injection to trigger the POST-synaptic spike
WARM_UP=200 # ms
DELTA_T=10 #ms

def _get_current_trace(freq,delta_t,t_stop,pre=False,test=False) :
    trace = np.zeros(t_stop/DT)
    for i in range(NO_REPS) :
        if(pre) :
            start_t = (0 + i* (1000.0/freq) + WARM_UP) 
        else :
            start_t = (0 + delta_t + i* (1000.0/freq) + WARM_UP) 
        end_t = (start_t+2)
        if(test) :
            print 'start_t=%g, end_t=%g (t_stop=%g, len(trace)=%f)' % (start_t,end_t,t_stop,len(trace))
        trace[start_t/DT:end_t/DT] = POST_AMP
    return trace

def test_pairing(delta_t=10,freqs=[0.1,10,20,40,50]) :
    gus = []
    for freq in freqs:#[0.1,10,20,40,50] :
        soma = h.Section()
        soma.insert('hh')
        soma.nseg = 1

        syn = h.STDPSynCC(soma(0.5))
        syn.delay_steps = 5.0 / DT # AP duration in raw number of integration steps
        syn.gbar=0.05

        # pre
        stim = h.NetStim()
        stim.number = NO_REPS
        stim.interval = 1000.0 / freq
        stim.start = WARM_UP # reverse: +DELTA_T
        stim.noise= 0
        nc = h.NetCon(stim,syn,0,0,0.1)

        #post
        ic = h.IClamp(soma(0.5))
        ic.delay = 0
        ic.dur=1e9
        total_time = 200+NO_REPS*(1000.0/freq)+100
        #print 'testing freq F=%g (t_stop=%i)' % (freq,total_time)
        current_trace = _get_current_trace(freq,delta_t,total_time,pre=False)
        current_vec = h.Vector(current_trace)
        current_vec.play(ic._ref_amp,DT)

        # set up some recording vectors
        trec = h.Vector()
        trec.record(h._ref_t)
        vrec = h.Vector()
        vrec.record(soma(0.5)._ref_v)
        grec = h.Vector()
        grec.record(syn._ref_g)
        gurec = h.Vector()
        gurec.record(syn._ref_g_update)
        gbrec = h.Vector()
        gbrec.record(syn._ref_gbar)
        um1s_rec = h.Vector()
        um1s_rec.record(syn._ref_um1s) 
        um2s_rec = h.Vector()
        um2s_rec.record(syn._ref_um2s) 
        #### END Of recording vectors

        h.dt = DT
        h.finitialize(-65)
        neuron.run(total_time)

        # not all arrays are required though...
        t = np.array(trec)
        v = np.array(vrec)
        g = np.array(grec)
        gu = np.array(gurec)
        gb = np.array(gbrec)
        um1s = np.array(um1s_rec)
        um2s = np.array(um2s_rec)

        # plt.figure()
        # plt.plot(t,v,'r',label='v')
        # plt.plot(t,um1s,'g',label='u_m1_s')
        # plt.plot(t,um2s,'b',label='u_m2_s')
        # plt.legend(loc=0)

        outcome=(gb[-1] - gb[0])*15.0/gb[0]*100.0
        print 'normalized weight change at F=%g : %f' % (freq,outcome)
        gus.append(outcome)

        # delete it all
        del(soma)
        del(syn)
        del(stim)
        del(ic)
        del(nc)
        del(trec);del(vrec);del(grec);del(gurec);del(gbrec)
    return gus

freqs =  [0.1,10,20,40,50]
pre_post_gus = test_pairing(delta_t=10,freqs=freqs)
post_pre_gus = test_pairing(delta_t=-10,freqs=freqs)

plt.figure(101)
plt.plot(freqs,pre_post_gus,label='pre-post')
plt.plot(freqs,post_pre_gus,label='post-pre')
plt.ylabel('normalized weight change(\%)')
plt.xlabel('rho (Hz)')
plt.suptitle('15 x 5  pairings at various rho='+str(freqs)+'Hz')
plt.legend(loc=0)
plt.show()
