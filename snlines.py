#!/usr/bin/env python
import optparse
import numpy as np
import pylab as py
import cPickle as pickle
from collections import deque

#########################################################
# Code for spectral line identification
# overplots a spectrum 
#########################################################

# some useful naming data
romannums = ('I','II','III', 'IV')
elements = ('H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni','Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'I', 'Te', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa', 'U')

py.ion()
pklfile = "kurucz.pkl"

# get command line arguments
parser = optparse.OptionParser()
parser.add_option("-v",dest="v_init")
parser.add_option("-t",dest="t_init")
parser.add_option("-n",dest="n_init")
parser.add_option("--xr",dest="xrange")
parser.add_option("--yr",dest="yrange")

# set defaults, if not specified
(opts, args) = parser.parse_args()
if (opts.v_init): vel  = float(opts.v_init)
else: vel = 1e4
if (opts.t_init): temp = float(opts.t_init)
else: temp = 1e4
if (opts.n_init): nshow = int(opts.n_init)
else: nshow = 1

# check for spectrum file
if (len(args)  == 0):
    print "You need to specify a spectrum file on the command line"
    exit()

# read input spectrum
xspec, yspec  = np.loadtxt(args[0], unpack = True)
yrange = (0,1.1*max(yspec))
xrange = (min(xspec),max(xspec))

# override the data ranges if set from command line
if (opts.xrange):
   xx = opts.xrange.split(',')
   x1 = float(xx[0])
   x2 = float(xx[1])
   xrange = [x1,x2]

if (opts.yrange):
   yy = opts.yrange.split(',')
   y1 = float(yy[0])
   y2 = float(yy[1])
   yrange = [y1,y2]

# read in line data, pickled by pickedata.py
with open(pklfile, 'rb') as f:
    linedict = pickle.load(f)

# plot it up
print "Using temperature " + str(temp) + " K",
print " and velocity " + str(vel) + " cm/s"
print "go for it (press ? for help, q to quit)"

# we start with no species
states = {}
specids = deque()

while (1):

    py.clf()
    py.plot(xspec,yspec)

    # make sure that things are set ok for active region
    try:
        active = specids[-1]
    except IndexError:
        active = None

    # overplot lines
    for id in specids:

        if id == active:
            color = 'red'
        else:
            color = 'black'

        state = states[id]

        for lam in state['lam'][:state['nshow']]:
            lam_obs = lam*(1 - state['vel']/3e5)
            x = [lam_obs,lam_obs]
            
            py.plot(x,yrange,color=color)

    py.xlabel('wavelength')
    py.ylabel('flux')
    if active is not None:
        py.title(states[active]['name'] + '; vel = ' + str(states[active]['vel']) + ' km/s')
             
    py.ylim(yrange)
    py.xlim(xrange)
    py.show()

    # get next command and break into command and args
    do = raw_input(">")
    try:
        docmd,doargs = do.split(' ', 1)
    except ValueError:
        docmd = do
        doargs = None

    if (docmd == 'q'): break;

    if (docmd == 'n'):
        e, i = doargs.split(' ')
        newid = (int(e),int(i))
        if newid in linedict:
            specids.append(newid)

            # start with one line showing 
            # and put it at the velocity of the previous active line
            if active is None:
                newvel = vel
            else:
                newvel = states[active]['vel']


            name = " ".join((elements[newid[0]-1], 
                             romannums[newid[1]]))

            
            # here are lines we even want to consider

            lam = linedict[newid]['lam']
            tau = linedict[newid]['tau']
            gf = linedict[newid]['gf']
            El = linedict[newid]['El']

            lamslice = (lam > xrange[0]) & (lam < xrange[1])

            # shift from reference temperature to temp
            k_B = 8.61733e-5
            tref = 1e4
            El = El * np.exp(El/k_B/tref) * np.exp(-El/k_B/temp)

            # pack this all in a dictionary
            states[newid] = {'vel': newvel, 
                             'nshow': 1,
                             'name': name,
                             'lam': lam[lamslice], 
                             'tau': tau[lamslice], 
                             'gf': gf[lamslice], 
                             'El' : El[lamslice]}
        else:
            print("No such line found in {}".format(pklfile))


    if (docmd == 'c'):
        try:
            nr = int(doargs)
        except (TypeError, ValueError):
            nr = 1
        specids.rotate(nr)

    if (docmd == 'k') and len(specids) >= 1:
        specids.pop()
                       
    if (docmd == 'd') and (active is not None): 
        states[active]['vel'] += 1e3
        print('velocity = ' + str(states[active]['vel']) + ' km/s')
    if (docmd == 'f') and (active is not None): 
        states[active]['vel'] -= 1e3
        print('velocity = ' + str(states[active]['vel']) + ' km/s')
    if (docmd == 'v') and (active is not None):
        states[active]['vel'] = float(input("new velocity (in km/s): "))
        
    if (docmd == 'a') and (active is not None): 

        try:
            nl = int(doargs)
        except (TypeError, ValueError):
            nl = 1

        states[active]['nshow'] = min(states[active]['nshow']+nl,
                                      len(states[active]['lam']))

    if (docmd == 'r') and (active is not None):

        try:
            nl = int(doargs)
        except (TypeError, ValueError):
            nl = 1

        states[active]['nshow'] = max(1,states[active]['nshow']-nl)

    if (docmd == 'e'):
        for id in specids:
            print("%4s %10.3e" % (states[id]['name'], states[id]['vel']))

    if (docmd == 'l'):
        print "   lambda        gf       E_low      tau"
        for id in specids:
            print("%4s @ vel = %10.3e km/s" % (states[id]['name'], states[id]['vel']))
            for i in range(states[id]['nshow']):
                print(" %10.3e %10.3e %10.3e %10.3e" % (states[id]['lam'][i],
                                                        states[id]['gf'][i],
                                                        states[id]['El'][i],
                                                        states[id]['tau'][i]))


    if (docmd == '?'):
        print("n #1 #2 = add element #1 with ionization state #2")
        print("c (#) = cycle the active species")
        print("k = remove active species")
        print("e = list all species")

        print("a (#) = add line(s) to active species")
        print("r (#) = remove line(s) from active species")
        print("l = list all current lines")

        print("d = increase velocity (blueshift)")
        print("f = decrease velocity (redshift)")
        print("v = reset the velocity by hand")

        print("q = quit")
        






