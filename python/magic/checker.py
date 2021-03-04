# -*- coding: utf-8 -*-
from magic import MagicTs, avgField, MagicSpectrum, MagicRadial, scanDir, MagicSetup
import numpy as np
from scipy.signal import argrelextrema
import matplotlib.pyplot as plt
import os

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[34m'
    OKGREEN = '\033[32m'
    MODERATE = '\033[33m'
    WARNING = '\033[91m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def MagicCheck(tstart=None):
    """
    This function is used to compute several sanity checks that can be evaluated
    if the power.TAG and some spectra have been produced in the current directory.
    If in addition the tInitAvg file is also there in the directory it averages
    only from this starting time.
    
    >>> MagicCheck(tstart=10.)
    """

    if os.path.exists('tInitAvg'):
        file = open('tInitAvg', 'r')
        tstart = float(file.readline().strip('\n'))
        file.close()

    if tstart is None:
        tstart = 0.

    ts = MagicTs(field='power', all=True, iplot=False)
    ts1 = MagicTs(field='dtE', all=True, iplot=False)
    mask = ( ts.time >= tstart )
    # Not super accurate if n_log_step changed but better that nothing
    n_steps = len(ts.time[mask]) * ts.n_log_step

    dEdt = ts.buoPower+ts.buoPower_chem+ts.ohmDiss+ts.viscDiss

    # Power balance:
    buoPower_avg = avgField(ts.time[mask], ts.buoPower[mask])
    buoPower_chem_avg = avgField(ts.time[mask], ts.buoPower_chem[mask])
    ohmDiss_avg = avgField(ts.time[mask], ts.ohmDiss[mask])
    viscDiss_avg = avgField(ts.time[mask], ts.viscDiss[mask])
    ratio = 100*abs(buoPower_avg+buoPower_chem_avg+ohmDiss_avg+viscDiss_avg) / \
            (buoPower_avg+buoPower_chem_avg) 

    print(bcolors.BOLD + bcolors.UNDERLINE + 'Power balance:' + bcolors.ENDC)
    print('Power injected   : {:.5e}'.format(buoPower_avg+buoPower_chem_avg))
    print('Power dissipated : {:.5e}'.format(-ohmDiss_avg-viscDiss_avg))
    st = 'Power mis-balance: {:.3f} %%'.format(ratio)
    if ratio <= 0.5:
        print(bcolors.OKGREEN + st + bcolors.ENDC)
    elif ratio > 0.5 and ratio <= 1.:
        print(bcolors.MODERATE + st + bcolors.ENDC)
    elif ratio > 1.:
        print(bcolors.WARNING + st + bcolors.ENDC)

    # Spikes catcher in the time series of power/dEdt
    print('\n' + bcolors.BOLD + bcolors.UNDERLINE + 'Time resolution:' + bcolors.ENDC)
    absdEdt_avg = avgField(ts1.time[mask], abs(ts1.dEdt[mask]))
    field = abs(dEdt[mask]-ts1.dEdt[mask])
    mask_spike = (field >= absdEdt_avg)
    l_spikes = False
    if mask_spike.any():
        l_spikes = True
    if l_spikes:
        print(bcolors.MODERATE + 'Sudden variations detected in power balance!' + bcolors.ENDC)
        ones = np.ones_like(ts.time[mask])
        ones[~mask_spike] = 0.
        ttot_spikes = np.trapz(ones, ts.time[mask])
        ttot = ts.time[-1]-ts.time[0]
        ratio = ttot_spikes/ttot
        print('    -Time fraction with spikes: {:.3f} %%'.format(100.*ratio))
        largest = abs(field).max()/absdEdt_avg
        st = '    -Largest event            : {:.2f} <|dE/dt|>'.format(largest)
        if largest > 10:
            print(bcolors.WARNING + st + bcolors.ENDC)
        else:
            print(bcolors.MODERATE + st + bcolors.ENDC)
    else:
        print(bcolors.OKGREEN + 'Power balance clean!' + bcolors.ENDC)

    # Timestep change occurence
    files = scanDir('timestep.*')
    if len(files) > 0:
        ts = MagicTs(field='timestep', iplot=False, all=True)
        mask = (ts.time >= tstart)
        ddt = np.diff(ts.time[mask])
        ddt_neq_zero = (ddt != 0.)
        ddt = ddt[ddt_neq_zero]
        print('\nNumber of time step changes   : {}'.format(len(ddt)))
        print('Number of iterations          : {}'.format(n_steps))
        freq = int( float(n_steps) / float(len(ddt)) )
        print('Average number of iterations with fixed time step size: {}'.format(freq))
        time = ts.time[mask][1:][ddt_neq_zero]
        dt = ts.dt[mask][1:][ddt_neq_zero]
        dtMean = avgField(time, dt)
        mask_changes = ( ddt <= 50*dtMean )
        ones = np.ones_like(time)
        ones[~mask_changes] = 0.
        ttot_changes = np.sum(ones*ddt)
        fast_change_ratio = 100*ttot_changes/(time[-1]-time[0])
        st = 'Fraction of time with frequent timestep changes (< 50 steps): {:.2f} %%'.format(
             fast_change_ratio)
        if fast_change_ratio < 2:
            print(bcolors.OKGREEN + st + bcolors.ENDC)
        elif fast_change_ratio >= 2 and fast_change_ratio <= 10:
            print(bcolors.MODERATE + st + bcolors.ENDC)
            print(bcolors.MODERATE + 'Maybe increase Courant factors!' + bcolors.ENDC)
        else:
            print(bcolors.WARNING + st + bcolors.ENDC)
            print(bcolors.WARNING + 'Probably increase Courant factors!' + bcolors.ENDC)

    # Dissipation lengthscales
    ts = MagicTs(field='par', all=True, iplot=False)
    mask = (ts.time >= tstart)
    lbDiss_avg = avgField(ts.time[mask], ts.lbDiss[mask])
    lvDiss_avg = avgField(ts.time[mask], ts.lvDiss[mask])

    ri = ts.radratio / (1.-ts.radratio)
    ro = 1. / (1.-ts.radratio)
    dmean = 0.5*(ri+ro)
    lTrunc = dmean*np.pi/ts.l_max

    print('\n' + bcolors.BOLD + bcolors.UNDERLINE + 'Angular resolution:' + bcolors.ENDC)
    print('Viscous dissipation scale: {:.3e}'.format(lvDiss_avg))
    if lbDiss_avg > 0:
        print('Ohmic dissipation scale  : {:.3e}'.format(lbDiss_avg))
    st = 'Angular truncation       : {:.3e}'.format(lTrunc)
    if lbDiss_avg > 0:
        lMin = min(lvDiss_avg, lbDiss_avg)
    else:
        lMin = lvDiss_avg
    ellMin = int(np.pi*dmean / lMin)
    nphi = int(3. * ellMin)
    if lTrunc < lMin:
        print(bcolors.OKGREEN + st + bcolors.ENDC)
    elif lTrunc >= lMin and lTrunc < 1.5*lMin:
        print(bcolors.MODERATE + st + bcolors.WARNING)
        st = 'You might need l_max={}, N_phi={}'.format(ellMin, nphi)
        print(bcolors.MODERATE + st + bcolors.WARNING)
    else:
        print(bcolors.WARNING + st + bcolors.WARNING)
        st = 'You might need l_max={}, N_phi={}'.format(ellMin, nphi)
        print(bcolors.WARNING + st + bcolors.WARNING)

    # Spectra
    dats = scanDir('kin_spec_ave.*')
    if len(dats) > 0:
        ave = True
    else:
        ave = False
    sp = MagicSpectrum(field='kin', iplot=False, ave=ave, quiet=True)
    ekin_l = sp.ekin_poll+sp.ekin_torl
    ratio = ekin_l.max()/ekin_l[-2]
    st = 'Vol. kin. energy spectra (largest/smallest): {:.2e}'.format(ratio)
    if ts.mode != 1:
        sp = MagicSpectrum(field='mag', iplot=False, ave=ave, quiet=True)
        emag_l = sp.emag_poll+sp.emag_torl
        ratio_mag = emag_l[2:].max()/emag_l[-2]
        ratio_cmb = sp.emagcmb_l[2:].max()/sp.emagcmb_l[-2]

        st_mag = 'Vol. mag. energy spectra (largest/smallest): {:.2e}'.format(ratio_mag)
        st_cmb = 'CMB mag. energy spectra (largest/smallest) : {:.2e}'.format(ratio_cmb)

    if ratio > 100:
        print(bcolors.OKGREEN + st + bcolors.ENDC)
    elif ratio <= 100 and ratio > 50:
        print(bcolors.MODERATE + st + bcolors.ENDC)
    else:
        print(bcolors.WARNING + st + bcolors.ENDC)
    if ts.mode != 1:
        if ratio_mag > 100:
            print(bcolors.OKGREEN + st_mag + bcolors.ENDC)
        elif ratio_mag <= 100 and ratio_mag > 50:
            print(bcolors.MODERATE + st_mag + bcolors.ENDC)
        else:
            print(bcolors.WARNING + st_mag + bcolors.ENDC)
        if ratio_cmb > 100:
            print(bcolors.OKGREEN + st_cmb + bcolors.ENDC)
        elif ratio_cmb <= 100 and ratio_cmb > 50:
            print(bcolors.MODERATE + st_cmb + bcolors.ENDC)
        else:
            print(bcolors.WARNING + st_cmb + bcolors.ENDC)

    # determine the relevant tags
    logs = scanDir('log.*')
    tags = []
    for lg in logs:
        stp = MagicSetup(nml=lg, quiet=True)
        if stp.start_time >= tstart and os.path.exists('eKinR.{}'.format(stp.tag)):
            tags.append(stp.tag)

    if len(tags) > 0:
        rad = MagicRadial(field='eKinR', iplot=False, quiet=True, tags=tags)
    else:
        rad = MagicRadial(field='eKinR', iplot=False, quiet=True)

    # Number of points in viscous BL
    print('\n' + bcolors.BOLD + bcolors.UNDERLINE + 'Radial resolution:' + bcolors.ENDC)
    if rad.ktopv != 1 and rad.kbotv != 1:
        eKR = rad.ekin_pol+rad.ekin_tor
    else:
        eKR = rad.ekin_pol
    ind = argrelextrema(eKR, np.greater)[0]
    ntop = ind[0]+1
    nbot = len(eKR)-ind[-1]
    nmin = min(nbot, ntop)
    st_bot = 'Number of points in bottom viscous B.L.: {}'.format(nbot)
    st_top = 'Number of points in top viscous B.L.   : {}'.format(ntop)

    if nbot >= 10:
        print(bcolors.OKGREEN + st_bot + bcolors.ENDC)
    elif nbot >=5 and nbot < 10:
        print(bcolors.MODERATE + st_bot + bcolors.ENDC)
    else:
        print(bcolors.WARNING + st_bot + bcolors.ENDC)
    if ntop >= 10:
        print(bcolors.OKGREEN + st_top + bcolors.ENDC)
    elif ntop >=5 and ntop < 10:
        print(bcolors.MODERATE + st_top + bcolors.ENDC)
    else:
        print(bcolors.WARNING + st_top + bcolors.ENDC)

    # Number of points in thermal BL
    if len(tags) > 0:
        rad = MagicRadial(field='heatR', iplot=False, quiet=True, tags=tags)
    else:
        rad = MagicRadial(field='heatR', iplot=False, quiet=True)
    if rad.ktops == 1:
        ind = argrelextrema(rad.entropy_SD, np.greater)[0]
        ntop = ind[0]+1
        st_top = 'Number of points in top thermal B.L.   : {}'.format(ntop)
        if ntop >= 10:
            print(bcolors.OKGREEN + st_top + bcolors.ENDC)
        elif ntop >=5 and ntop < 10:
            print(bcolors.MODERATE + st_top + bcolors.ENDC)
        else:
            print(bcolors.WARNING + st_top + bcolors.ENDC)

    if rad.kbots == 1:
        ind = argrelextrema(rad.entropy_SD, np.greater)[0]
        nbot = len(rad.radius)-ind[-1]
        st_bot = 'Number of points in bottom thermal B.L.: {}'.format(nbot)
        if nbot >= 10:
            print(bcolors.OKGREEN + st_bot + bcolors.ENDC)
        elif nbot >=5 and nbot < 10:
            print(bcolors.MODERATE + st_bot + bcolors.ENDC)
        else:
            print(bcolors.WARNING + st_bot + bcolors.ENDC)
