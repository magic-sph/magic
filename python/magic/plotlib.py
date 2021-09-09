# -*- coding: utf-8 -*-
from magic.setup import defaultLevels, defaultCm
import numpy as np
import matplotlib.pyplot as plt


def hammer2cart(ttheta, pphi, colat=False):
    """
    This function is used to define the Hammer projection used when
    plotting surface contours in :py:class:`magic.Surf`

    >>> # Load Graphic file
    >>> gr = MagicGraph()
    >>> # Meshgrid
    >>> pphi, ttheta = mgrid[-np.pi:np.pi:gr.nphi*1j, np.pi/2.:-np.pi/2.:gr.ntheta*1j]
    >>> x,y = hammer2cart(ttheta, pphi)
    >>> # Contour plots
    >>> contourf(x, y, gr.vphi)

    :param ttheta: meshgrid [nphi, ntheta] for the latitudinal direction
    :type ttheta: numpy.ndarray
    :param pphi: meshgrid [nphi, ntheta] for the azimuthal direction
    :param colat: colatitudes (when not specified a regular grid is assumed)
    :type colat: numpy.ndarray
    :returns: a tuple that contains two [nphi, ntheta] arrays: the x, y meshgrid
              used in contour plots
    :rtype: (numpy.ndarray, numpy.ndarray)
    """

    if not colat: # for lat and phi \in [-pi, pi]
        xx = 2.*np.sqrt(2.) * np.cos(ttheta)*np.sin(pphi/2.)\
             /np.sqrt(1.+np.cos(ttheta)*np.cos(pphi/2.))
        yy = np.sqrt(2.) * np.sin(ttheta)\
             /np.sqrt(1.+np.cos(ttheta)*np.cos(pphi/2.))
    else:  # for colat and phi \in [0, 2pi]
        xx = -2.*np.sqrt(2.) * np.sin(ttheta)*np.cos(pphi/2.)\
             /np.sqrt(1.+np.sin(ttheta)*np.sin(pphi/2.))
        yy = np.sqrt(2.) * np.cos(ttheta)\
             /np.sqrt(1.+np.sin(ttheta)*np.sin(pphi/2.))
    return xx, yy


def cut(dat, vmax=None, vmin=None):
    """
    This functions truncates the values of an input array that are beyond
    vmax or below vmin and replace them by vmax and vmin, respectively.

    >>> # Keep only values between -1e3 and 1e3
    >>> datNew = cut(dat, vmin=-1e3, vmax=1e3)

    :param dat: an input array
    :type dat: numpy.ndarray
    :param vmax: maximum upper bound
    :type vmax: float
    :param vmin: minimum lower bound
    :type vmin: float
    :returns: an array where the values >=vmax have been replaced by vmax
              and the values <=vmin have been replaced by vmin
    :rtype: numpy.ndarray
    """
    if vmax is not None:
        mask = np.where(dat>=vmax, 1, 0)
        dat = dat*(mask == 0) + vmax*(mask == 1)
        normed = False
    if vmin is not None:
        mask = np.where(dat<=vmin, 1, 0)
        dat = dat*(mask == 0) + vmin*(mask == 1)
        normed = False
    return dat


def equatContour(data, radius, minc=1, label=None, levels=defaultLevels,
                 cm=defaultCm, normed=True, vmax=None, vmin=None, cbar=True,
                 tit=True, normRad=False, deminc=True, bounds=True):
    """
    Plot the equatorial cut of a given field

    :param data: the input data (an array of size (nphi,nr))
    :type data: numpy.ndarray
    :param radius: the input radius
    :type radius: numpy.ndarray
    :param minc: azimuthal symmetry
    :type minc: int
    :param label: the name of the input physical quantity you want to
                  display
    :type label: str
    :param normRad: when set to True, the contour levels are normalised
                    radius by radius (default is False)
    :type normRad: bool
    :param levels: the number of levels in the contour
    :type levels: int
    :param cm: name of the colormap ('jet', 'seismic', 'RdYlBu_r', etc.)
    :type cm: str
    :param tit: display the title of the figure when set to True
    :type tit: bool
    :param cbar: display the colorbar when set to True
    :type cbar: bool
    :param vmax: maximum value of the contour levels
    :type vmax: float
    :param vmin: minimum value of the contour levels
    :type vmin: float
    :param normed: when set to True, the colormap is centered around zero.
                   Default is True, except for entropy/temperature plots.
    :type normed: bool
    :param deminc: a logical to indicate if one wants do get rid of the
                   possible azimuthal symmetry
    :type deminc: bool
    :param bounds: a boolean to determine if one wants to plot the limits
                   of the domain (True by default)
    :type bounds: bool
    """

    nphi, ntheta = data.shape

    if deminc:
        phi = np.linspace(0., 2.*np.pi, nphi)
    else:
        phi = np.linspace(0., 2.*np.pi/minc, nphi)
    rr, pphi = np.meshgrid(radius, phi)
    xx = rr * np.cos(pphi)
    yy = rr * np.sin(pphi)

    if normRad: # Normalise each radius
        maxS = np.sqrt(np.mean(data**2, axis=0))
        data[:, maxS!=0.] /= maxS[maxS!=0.]

    if tit and label is not None:
        if cbar:
            fig = plt.figure(figsize=(6.5,5.5))
            ax = fig.add_axes([0.01, 0.01, 0.76, 0.9])
        else:
            fig = plt.figure(figsize=(5,5.5))
            ax = fig.add_axes([0.01, 0.01, 0.98, 0.9])
        ax.set_title(label, fontsize=24)
    else:
        if cbar:
            fig = plt.figure(figsize=(6.5,5))
            ax = fig.add_axes([0.01, 0.01, 0.76, 0.98])
        else:
            fig = plt.figure(figsize=(5, 5))
            ax = fig.add_axes([0.01, 0.01, 0.98, 0.98])

    cmap = plt.get_cmap(cm)
    if vmax is not None or vmin is not None:
        normed = False
        cs = np.linspace(vmin, vmax, levels)
        im = ax.contourf(xx, yy, data, cs, cmap=cmap, extend='both')
    else:
        cs = levels
        im = ax.contourf(xx, yy, data, cs, cmap=cmap)
        #im = ax.pcolormesh(xx, yy, data, cmap=cmap, antialiased=True)

    if bounds:
        ax.plot(radius[0]*np.cos(phi), radius[0]*np.sin(phi), 'k-', lw=1.5)
        ax.plot(radius[-1]*np.cos(phi), radius[-1]*np.sin(phi), 'k-', lw=1.5)

        if not deminc and minc > 1:
            ax.plot(radius, np.zeros_like(radius), 'k-', lw=1.5)
            xa = radius[-1]*np.cos(2.*np.pi/minc)
            ya = radius[-1]*np.sin(2.*np.pi/minc)
            xb = radius[0]*np.cos(2.*np.pi/minc)
            x = np.linspace(xa, xb, 32)
            y = np.tan(2.*np.pi/minc)*(x-xa)+ya
            ax.plot(x, y, 'k-', lw=1.5)
            ax.plot(radius, np.zeros_like(radius), 'k-', lw=1.5)

    eps = 1e-4
    if xx.min() < -eps:
        ax.set_xlim(1.01*xx.min(), 1.01*xx.max())
    elif abs(xx.min()) < eps :
        ax.set_xlim(xx.min()-0.01, 1.01*xx.max())
    else:
        ax.set_xlim(0.99*xx.min(), 1.01*xx.max())
    if yy.min() < -eps:
        ax.set_ylim(1.01*yy.min(), 1.01*yy.max())
    elif abs(yy.min()) < eps:
        ax.set_ylim(yy.min()-0.01, 1.01*yy.max())
    else:
        ax.set_ylim(0.99*yy.min(), 1.01*yy.max())
    ax.axis('off')

    # Add the colorbar at the right place
    pos = ax.get_position()
    l, b, w, h = pos.bounds
    if cbar:
        if tit and label is not None:
            cax = fig.add_axes([0.85, 0.46-0.7*h/2., 0.03, 0.7*h])
        else:
            cax = fig.add_axes([0.85, 0.5-0.7*h/2., 0.03, 0.7*h])
        mir = fig.colorbar(im, cax=cax)

    # Normalise data
    if normed:
        im.set_clim(-max(abs(data.max()), abs(data.min())),
                     max(abs(data.max()), abs(data.min())))

    #To avoid white lines on pdfs

    for c in im.collections:
        c.set_edgecolor("face")

    return fig, xx, yy


def merContour(data, radius, label=None, levels=defaultLevels, cm=defaultCm,
               normed=True, vmax=None, vmin=None, cbar=True, tit=True,
               fig=None, ax=None, bounds=True, lines=False):
    """
    Plot a meridional cut of a given field

    :param data: the input data (an array of size (ntheta,nr))
    :type data: numpy.ndarray
    :param radius: the input radius
    :type radius: numpy.ndarray
    :param label: the name of the input physical quantity you want to
                  display
    :type label: str
    :param levels: the number of levels in the contour
    :type levels: int
    :param cm: name of the colormap ('jet', 'seismic', 'RdYlBu_r', etc.)
    :type cm: str
    :param tit: display the title of the figure when set to True
    :type tit: bool
    :param cbar: display the colorbar when set to True
    :type cbar: bool
    :param vmax: maximum value of the contour levels
    :type vmax: float
    :param vmin: minimum value of the contour levels
    :type vmin: float
    :param normed: when set to True, the colormap is centered around zero.
                   Default is True, except for entropy/temperature plots.
    :type normed: bool
    :param bounds: a boolean to determine if one wants to plot the limits
                   of the domain (True by default)
    :type bounds: bool
    :param fig: a pre-existing figure (if needed)
    :type fig: matplotlib.figure.Figure
    :param ax: a pre-existing axis
    :type ax: matplotlib.axes._subplots.AxesSubplot
    :param lines: when set to True, over-plot solid lines to highlight
                  the limits between two adjacent contour levels
    :type lines: bool
    """
    ntheta, nr = data.shape

    th = np.linspace(0, np.pi, ntheta)
    rr, tth = np.meshgrid(radius, th)
    xx = rr * np.sin(tth)
    yy = rr * np.cos(tth)

    if tit and label is not None:
        if cbar:
            fsize = (5, 7.5)
            bb = [0.01, 0.01, 0.69, 0.91]
        else:
            fsize = (3.5, 7.5)
            bb = [0.01, 0.01, 0.98, 0.91]
    else:
        if cbar:
            fsize = (5, 7)
            bb = [0.01, 0.01, 0.69, 0.98]
        else:
            fsize = (3.5, 7)
            bb = [0.01, 0.01, 0.98, 0.98]

    if fig is None:
        fig = plt.figure(figsize=fsize)
        if ax is None:
            ax = fig.add_axes(bb)

    if tit and label is not None:
        ax.set_title(label, fontsize=24)

    cmap = plt.get_cmap(cm)
    if vmax is not None and vmin is not None:
        normed = False
        cs = np.linspace(vmin, vmax, levels)
        im = ax.contourf(xx, yy, data, cs, cmap=cmap, extend='both')
        if lines:
            ax.contour(xx, yy, data, cs, colors=['k'], linewidths=0.5,
                       extend='both', linestyles=['-'])
    else:
        cs = levels
        im = ax.contourf(xx, yy, data, cs, cmap=cmap)
        if lines:
            ax.contour(xx, yy, data, cs, colors=['k'], linewidths=0.5,
                       linestyles=['-'])
        #im = ax.pcolormesh(xx, yy, data, cmap=cmap, antialiased=True)

    #To avoid white lines on pdfs

    for c in im.collections:
        c.set_edgecolor("face")

    if bounds:
        ax.plot((radius[0])*np.sin(th), (radius[0])*np.cos(th), 'k-')
        ax.plot((radius[-1])*np.sin(th), (radius[-1])*np.cos(th), 'k-')
        ax.plot([0., 0.], [radius[-1], radius[0]], 'k-')
        ax.plot([0., 0.], [-radius[-1], -radius[0]], 'k-')

    eps = 1e-4
    if xx.min() < -eps:
        ax.set_xlim(1.01*xx.min(), 1.01*xx.max())
    elif abs(xx.min()) < eps :
        ax.set_xlim(xx.min()-0.01, 1.01*xx.max())
    else:
        ax.set_xlim(0.99*xx.min(), 1.01*xx.max())
    if yy.min() < -eps:
        ax.set_ylim(1.01*yy.min(), 1.01*yy.max())
    elif abs(yy.min()) < eps:
        ax.set_ylim(yy.min()-0.01, 1.01*yy.max())
    else:
        ax.set_ylim(0.99*yy.min(), 1.01*yy.max())
    ax.axis('off')

    # Add the colorbar at the right place
    pos = ax.get_position()
    l, b, w, h = pos.bounds
    if cbar:
        if tit and label is not None:
            cax = fig.add_axes([0.82, 0.46-0.7*h/2., 0.04, 0.7*h])
        else:
            cax = fig.add_axes([0.82, 0.5-0.7*h/2., 0.04, 0.7*h])
        mir = fig.colorbar(im, cax=cax)

    # Normalisation of the contours around zero
    if normed is True:
        im.set_clim(-max(abs(data.max()), abs(data.min())),
                     max(abs(data.max()), abs(data.min())))

    return fig, xx, yy, im


def radialContour(data, rad=0.85, label=None, proj='hammer', lon_0=0., vmax=None,
                  vmin=None, lat_0=30., levels=defaultLevels, cm=defaultCm,
                  normed=True, cbar=True, tit=True, lines=False, fig=None,
                  ax=None):
    """
    Plot the radial cut of a given field

    :param data: the input data (an array of size (nphi,ntheta))
    :type data: numpy.ndarray
    :param rad: the value of the selected radius
    :type rad: float
    :param label: the name of the input physical quantity you want to
                  display
    :type label: str
    :param proj: the type of projection. Default is Hammer, in case
                 you want to use 'ortho' or 'moll', then Basemap is
                 required.
    :type proj: str
    :param levels: the number of levels in the contour
    :type levels: int
    :param cm: name of the colormap ('jet', 'seismic', 'RdYlBu_r', etc.)
    :type cm: str
    :param tit: display the title of the figure when set to True
    :type tit: bool
    :param cbar: display the colorbar when set to True
    :type cbar: bool
    :param lines: when set to True, over-plot solid lines to highlight
                  the limits between two adjacent contour levels
    :type lines: bool
    :param vmax: maximum value of the contour levels
    :type vmax: float
    :param vmin: minimum value of the contour levels
    :type vmin: float
    :param normed: when set to True, the colormap is centered around zero.
                   Default is True, except for entropy/temperature plots.
    :type normed: bool
    :param fig: a pre-existing figure (if needed)
    :type fig: matplotlib.figure.Figure
    :param ax: a pre-existing axis
    :type ax: matplotlib.axes._subplots.AxesSubplot
    """

    nphi, ntheta = data.shape

    phi = np.linspace(-np.pi, np.pi, nphi)
    theta = np.linspace(np.pi/2, -np.pi/2, ntheta)
    pphi, ttheta = np.mgrid[-np.pi:np.pi:nphi*1j, np.pi/2.:-np.pi/2.:ntheta*1j]
    lon2 = pphi * 180./np.pi
    lat2 = ttheta * 180./np.pi

    circles = np.r_[-60., -30., 0., 30., 60.]
    delon = 60.
    meridians = np.arange(-180+delon, 180, delon)

    if fig is None and ax is None:
        if proj == 'moll' or proj == 'hammer':
            if tit and label is not None:
                if cbar:
                    fig = plt.figure(figsize=(9,4.5))
                    ax = fig.add_axes([0.01, 0.01, 0.87, 0.87])
                else:
                    fig = plt.figure(figsize=(8,4.5))
                    ax = fig.add_axes([0.01, 0.01, 0.98, 0.87])
                ax.set_title('{}: r/ro = {:.3f}'.format(label, rad),
                             fontsize=24)
            else:
                if cbar:
                    fig = plt.figure(figsize=(9,4))
                    ax = fig.add_axes([0.01, 0.01, 0.87, 0.98])
                else:
                    fig = plt.figure(figsize=(8,4))
                    ax = fig.add_axes([0.01, 0.01, 0.98, 0.98])
                #tit1 = r'{:.2f} Ro'.format(rad)
                #ax.text(0.12, 0.9, tit1, fontsize=16,
                      #horizontalalignment='right',
                      #verticalalignment='center',
                      #transform = ax.transAxes)
        else:
            if tit and label is not None:
                if cbar:
                    fig = plt.figure(figsize=(6,5.5))
                    ax = fig.add_axes([0.01, 0.01, 0.82, 0.9])
                else:
                    fig = plt.figure(figsize=(5,5.5))
                    ax = fig.add_axes([0.01, 0.01, 0.98, 0.9])
                ax.set_title('{}: r/ro = {:.3f}'.format(label, rad),
                             fontsize=24)
            else:
                if cbar:
                    fig = plt.figure(figsize=(6,5))
                    ax = fig.add_axes([0.01, 0.01, 0.82, 0.98])
                else:
                    fig = plt.figure(figsize=(5,5))
                    ax = fig.add_axes([0.01, 0.01, 0.98, 0.98])
                tit1 = r'{:.2f} Ro'.format(rad)
                ax.text(0.12, 0.9, tit1, fontsize=16,
                      horizontalalignment='right',
                      verticalalignment='center',
                      transform = ax.transAxes)

    if proj != 'hammer':
        from mpl_toolkits.basemap import Basemap
        map = Basemap(projection=proj, lon_0=lon_0, lat_0=lat_0,
                      resolution='c')
        map.drawparallels([0.], dashes=[2, 3], linewidth=0.5)
        map.drawparallels(circles, dashes=[2,3], linewidth=0.5)
        map.drawmeridians(meridians, dashes=[2,3], linewidth=0.5)
        map.drawmeridians([-180], dashes=[20,0], linewidth=0.5)
        map.drawmeridians([180], dashes=[20,0], linewidth=0.5)
        x, y = list(map(lon2, lat2))
    else:
        x, y = hammer2cart(ttheta, pphi)
        for lat0 in circles:
            x0, y0 = hammer2cart(lat0*np.pi/180., phi)
            ax.plot(x0, y0, 'k:', linewidth=0.7)
        for lon0 in meridians:
            x0, y0 = hammer2cart(theta, lon0*np.pi/180.)
            ax.plot(x0, y0, 'k:', linewidth=0.7)
        xxout, yyout  = hammer2cart(theta, -np.pi-1e-3)
        xxin, yyin  = hammer2cart(theta, np.pi+1e-3)
        ax.plot(xxin, yyin, 'k-')
        ax.plot(xxout, yyout, 'k-')
        ax.axis('off')

    cmap = plt.get_cmap(cm)

    if proj == 'ortho':
        lats = np.linspace(-90., 90., ntheta)
        dat = map.transform_scalar(data.T, phi*180/np.pi, lats,
                                   nphi, ntheta, masked=True)
        im = map.imshow(dat, cmap=cmap)
    else:
        if vmax is not None or vmin is not None:
            normed = False
            cs = np.linspace(vmin, vmax, levels)
            im = ax.contourf(x, y, data, cs, cmap=cmap, extend='both')
            if lines:
                ax.contour(x, y, data, cs, colors=['k'], linewidths=0.5,
                           extend='both', linestyles=['-'])
                #ax.contour(x, y, data, 1, colors=['k'])
            #im = ax.pcolormesh(x, y, data, cmap=cmap, antialiased=True)
        else:
            cs = levels
            im = ax.contourf(x, y, data, cs, cmap=cmap)
            if lines:
                ax.contour(x, y, data, cs, colors=['k'], linewidths=0.5,
                           linestyles=['-'])
            #im = ax.pcolormesh(x, y, data, cmap=cmap, antialiased=True)


    # Add the colorbar at the right place
    pos = ax.get_position()
    l, b, w, h = pos.bounds
    if cbar:
        if tit and label is not None:
            cax = fig.add_axes([0.9, 0.46-0.7*h/2., 0.03, 0.7*h])
        else:
            cax = fig.add_axes([0.9, 0.51-0.7*h/2., 0.03, 0.7*h])
        mir = fig.colorbar(im, cax=cax)

    # Normalise around zero
    if normed:
        im.set_clim(-max(abs(data.max()), abs(data.min())),
                     max(abs(data.max()), abs(data.min())))

    #To avoid white lines on pdfs

    for c in im.collections:
        c.set_edgecolor("face")

    return fig
