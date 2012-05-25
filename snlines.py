#!/usr/bin/env python
import optparse
import numpy as np
import pylab as py

#########################################################
# Code for spectral line identification
# overplots a spectrum 
#########################################################


py.ion()
linefile = "kurucz_cd23_cut.dat"

# get command line arguments
parser = optparse.OptionParser()
parser.add_option("-v",dest="v_init")
parser.add_option("-t",dest="t_init")
parser.add_option("-n",dest="n_init")
parser.add_option("-z",dest="Z_init")
parser.add_option("-Z",dest="Z_init")
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

# Get element Z and ionization stage (0 is neutral)
if (opts.Z_init):
    zz = opts.Z_init.split('.')
    elem = int(zz[0])
    ion  = int(zz[1])
    print 'elem = ' + str(elem) + '; ion = ' + str(ion)
else:
    elem = int(input("element (atomic number Z): "))
    ion  = int(input("ion stage (0 for neutral): "))    

# read input spectrum
data  = np.loadtxt(args[0])
xspec = data[:,0]
yspec = data[:,1]
yrange = (0,1.1*max(yspec))
xrange = (min(xspec),max(xspec))

if (opts.xrange):
   xx = opts.xrange.split(',')
   x1 = float(xx[0])
   x2 = float(xx[1])
   xrange = [x1,x2]

if (opts.yrange):
   xx = opts.yrange.split(',')
   x1 = float(xx[0])
   x2 = float(xx[1])
   yrange = [y1,y2]


# read line list
f = open(linefile)
lam  = np.zeros(0)
gf   = np.zeros(0)
El   = np.zeros(0)
for line in f:
    data = line.split()
    e = int(data[3])
    i = int(data[4])
    if (e == elem and i == ion):
        l = float(data[0])
        if (l > xrange[0] and l < xrange[1]):
            lam = np.append(lam,l)
            gf  = np.append(gf,float(data[1]))
            El  = np.append(El,float(data[2]))

# check for line number
if (len(lam) == 0):
    print "No lines for that element found in this wavelength range"
    exit()

# compute (relative) sobolev optical depths assuming LTE
# (depends on oscilator strength and boltzmann factor)
k_B = 8.61733e-5
tau = lam*gf*np.exp(-El/k_B/temp)
tau = tau/tau.max()

# sort lines by sobolev optical depth
ind = tau.argsort()
ind = ind[::-1]
tau = tau[ind]
lam = lam[ind]
gf  = gf[ind]
El  = El[ind]


# plot it up
print "Using temperature " + str(temp) + " K",
print " and velocity " + str(vel) + " cm/s"
print "go for it (press ? for help, q to quit)"
while (1):

    py.clf()
    py.plot(xspec,yspec)

    # overplot lines
    for i in range(nshow):
        lam_obs = lam[i]*(1 - vel/3e5)
        x = [lam_obs,lam_obs]
        py.plot(x,yrange,color='black')

    py.xlabel('wavelength')
    py.ylabel('flux')
    py.title('vel = ' + str(vel) + ' km/s; Z = ' + str(elem) + '; ion = ' + str(ion))
             

    py.ylim(yrange)
    py.xlim(xrange)
    py.show()

    # get next command
    do = raw_input(">")

    if (do == 'q'): break;
    if (do == 'd'): 
        vel = vel + 1e3
        print 'velocity = ' + str(vel) + ' km/s'
    if (do == 'f'): 
        vel = vel - 1e3
        print 'velocity = ' + str(vel) + ' km/s'
    if (do == 'v'):
        vel = float(input("new velocity (in km/s): "))
        
    if (do == 'a'): 
        print "adding: %10.3e %10.3e %10.3e %10.3e" % (lam[i],gf[i],El[i],tau[i])
        nshow = nshow + 1
        
    if (do == 'r'): nshow = nshow - 1
    
    if (do == 'l'):
        print "   lambda        gf       E_low      tau"
        for i in range(nshow):
            print "%10.3e %10.3e %10.3e %10.3e" % (lam[i],gf[i],El[i],tau[i])

    if (do == '?'):
        print "a = add a line"
        print "r = remove a line"
        print "l = list all current lines"
        print "d = increase velocity (blueshift)"
        print "f = decrease velocity (blueshift)"
        print "v = reset the velocity by hand"
        print "q = quit"
        






