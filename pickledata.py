#!/usr/bin/env python
import optparse
import numpy as np
import cPickle as pickle

def main(linefile,pklfile):

    # read line list
    linedata = np.loadtxt(linefile, 
                          dtype = {'names' : ('lam','gf','El','e','i'), 
                                   'formats' : ('f4','f4','f4','i4','i4')})
    
    # here are all the species
    allelts = set([(e,i) for e,i in zip(linedata['e'],linedata['i'])])

    # go through and fill up the linedict
    linedict = {}
    for elt in allelts:

        # get an array slice for just those elements with this (e,i) pair
        slice = (linedata['e']==elt[0]) & (linedata['i']==elt[1])
        
        lam = linedata['lam'][slice]
        gf = linedata['gf'][slice]
        El = linedata['El'][slice]
        
        # compute (relative) sobolev optical depths assuming LTE
        # (depends on oscilator strength and boltzmann factor)
        k_B = 8.61733e-5
        tref = 1e4
        tau = lam*gf*np.exp(-El/k_B/tref)
        tau = tau/tau.max()
        
        # sort lines by sobolev optical depth
        ind = tau.argsort()
        ind = ind[::-1]
        tau = tau[ind]
        lam = lam[ind]
        gf  = gf[ind]
        El  = El[ind]
    
        # dump the physical data into a dictionary indexed by (e,i) tuple
        linedict[elt] = {'lam': lam, 'tau': tau, 'gf': gf, 'El' : El}

    # pickle our hard work
    with open(pklfile, 'wb') as f:
        pickle.dump(linedict, f)


if __name__ == "__main__":

    linefile = "kurucz_cd23_cut.dat"
    pklfile = "kurucz.pkl"

    main(linefile,pklfile)
