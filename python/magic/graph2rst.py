# -*- coding: utf-8 -*-
from . import MagicGraph, SpectralTransforms, thetaderavg, phideravg, \
              MagicRadial
import numpy as np
import os
import struct
import sys


def bytes_per_record(record_type):
    """
    This routine determines the number of bytes that corresponds
    to an input record

    :param record_type: the type of record
    :type record_type: str
    """

    if record_type == '>f' or record_type == 'f':
        return 4
    if record_type == '>d' or record_type == 'd':
        return 8
    if record_type == '>i' or record_type == 'i':
        return 4

    sys.stderr.write('unknown record_type\n')
    sys.exit(-1)


def write_record(f, record_type, array):
    """
    This subroutine writes an array into a binary file including
    the corresponding record markers to make it readable by Fortran

    :param f: the binary file
    :type f: file
    :param record_type: the array type (can be 'f', 'd', 'i')
    :type record_type: str
    :param array: the array to be written
    :type array: np.ndarray
    """

    # starting marker
    f.write(struct.pack('>i', bytes_per_record(record_type)*len(array)))

    for b in array:
        f.write(struct.pack(record_type, b))

    # read end marker
    f.write(struct.pack('>i', bytes_per_record(record_type)*len(array)))


def write_single_record(f, record_type, val):
    """
    This subroutine writes one variable into a binary file including
    the corresponding record markers to make it readable by Fortran

    :param f: the binary file
    :type f: file
    :param record_type: the type of the variable (can be 'f', 'd', 'i')
    :type record_type: str
    :param val: the scalar
    :type array: float or int
    """

    f.write(struct.pack('>i', bytes_per_record(record_type)))
    f.write(struct.pack(record_type, val))
    f.write(struct.pack('>i', bytes_per_record(record_type)))


class Graph2Rst:
    """
    This class allows to transform an input graphic file into a checkpoint
    file format that can be read by MagIC to restart a simulation.

    >>> # Load a Graphic File
    >>> gr = MagicGraph()
    >>> # Produce the file checkpoint_ave.from_G with dt=1e-4
    >>> Graph2Rst(gr, filename='checkpoint_ave.from_G', dt=1e-4):q

    .. note:: So far the endianness is hard-coded an set to "big_endian".
              This might become useful to have an extra option to control
              it as well.
    """

    def __init__(self, gr, filename='checkpoint_ave', dt=-1):
        """
        :param gr: the input graphic file one wants to convert into a restart
                   file
        :type gr: magic.MagicGraph
        :param filename: name of the checkpoint file
        :type filename: str
        :param dt: time step size
        :type dt: float
        """

        if hasattr(gr, 'tag'):
            tag = gr.tag

            if os.path.exists('anel.%s' % tag):
                r = MagicRadial(field='anel', iplot=False)
                rho0 = r.rho0
            else:
                rho0 = np.ones_like(gr.radius)

        else:
            rho0 = np.ones_like(gr.radius)

        f = open(filename, 'wb')

        version = 1

        write_single_record(f, '>i', version)

        time = gr.time.astype(np.float64)
        if dt < 0:
            dt = gr.dtMax

        f.write(struct.pack('>i', 8+8+4))  # 2 doubles + 1 int
        f.write(struct.pack('>d', time))
        f.write(struct.pack('>d', dt))
        f.write(struct.pack('>i', gr.n_time_steps))
        f.write(struct.pack('>i', 8+8+4))  # 2 doubles + 1 int

        ra = gr.ra.reshape((1))[0]
        pr = gr.pr.reshape((1))[0]
        if type(gr.raxi) == np.ndarray:
            raxi = gr.raxi.reshape((1))[0]
            sc = gr.sc.reshape((1))[0]
        else:
            raxi = gr.raxi
            sc = gr.sc
        prmag = gr.prmag.reshape((1))[0]
        ek = gr.ek.reshape((1))[0]
        radratio = gr.radratio.reshape((1))[0]
        par = np.r_[ra, pr, raxi, sc, prmag, ek, radratio, gr.sigma_ratio]
        par = par.astype(np.float64)
        write_record(f, '>d', par)

        # Grid resolution
        res = np.r_[gr.n_r_max, gr.n_theta_max, gr.n_phi_tot, gr.minc,
                    gr.nalias, gr.n_r_ic_max]
        write_record(f, '>i', res)

        # Radial scheme
        if gr.radial_scheme == 'CHEB':
            version = 'cheb'+'%68s' % ''
            n_max = gr.n_cheb_max
            if gr.l_newmap == 'F':
                o_bound = 0
            else:
                o_bound = 1
            alph1 = gr.alph1
            alph2 = gr.alph2
        else:
            version = 'fd'+'%70s' % ''
            n_max = gr.fd_order
            o_bound = gr.fd_order_bound
            alph1 = gr.fd_stretch
            alph2 = gr.fd_ratio

        s = bytes(version)
        f.write(struct.pack('>i', 8+8+4+4+len(s)))  # 2 doubles + 2 int + 1 str
        f.write(struct.pack(">%ds" % (len(s),), s))
        f.write(struct.pack('>i', n_max))
        f.write(struct.pack('>i', o_bound))
        f.write(struct.pack('>d', alph1))
        f.write(struct.pack('>d', alph2))
        f.write(struct.pack('>i', 8+8+4+4+len(s)))  # 2 doubles + 2 int + 1str

        # Torques: dummy
        dumm = np.zeros(15, dtype=np.float64)
        dumm[-1] = dt
        write_record(f, '>d', dumm)

        # Flags
        if gr.mode in [2, 3, 7, 8, 9, 10] or ra == 0.:
            l_heat = 0
        else:
            l_heat = 1

        if raxi > 0. or raxi < 0.:
            l_chem = 1
        else:
            l_chem = 0

        if gr.mode in [0, 2, 3, 6, 8, 9]:
            l_mag = 1
        else:
            l_mag = 0

        if gr.sigma_ratio == 0.:
            l_cond_ic = 0
        else:
            l_cond_ic = 1

        flags = np.r_[l_heat, l_chem, l_mag, l_cond_ic]
        write_record(f, '>i', flags)

        # At this stage we need to initiate the transforms

        sh = SpectralTransforms(l_max=gr.l_max, lm_max=gr.lm_max, minc=gr.minc,
                                n_theta_max=gr.n_theta_max)

        # Calculate and store the poloidal potential using vr
        pol = np.zeros((gr.n_r_max, gr.lm_max), dtype=np.complex128)
        for i in range(gr.n_r_max):
            vr = sh.spat_spec(gr.vr[:, :, i])
            pol[i, 1:] = vr[1:]/(sh.ell[1:]*(sh.ell[1:]+1)) * \
                         gr.radius[i]**2 * rho0[i]

        f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_max))
        pol.byteswap(True)
        pol.tofile(f)
        f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_max))

        pol = np.zeros((gr.n_r_max, gr.lm_max), dtype=np.complex128)
        f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_max))
        pol.byteswap(True)
        pol.tofile(f)
        f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_max))

        # Calculate the toroidal potential using wr
        tor = np.zeros((gr.n_r_max, gr.lm_max), dtype=np.complex128)

        th3D = np.zeros_like(gr.vr)
        rr3D = np.zeros_like(th3D)
        for i in range(gr.n_theta_max):
            th3D[:, i, :] = gr.colatitude[i]
        for i in range(gr.n_r_max):
            rr3D[:, :, i] = gr.radius[i]
        s3D = rr3D*np.sin(th3D)
        omr = 1./s3D*(thetaderavg(np.sin(th3D)*gr.vphi, order=4) -
                      phideravg(gr.vtheta, minc=gr.minc))

        for i in range(gr.n_r_max):
            om = sh.spat_spec(omr[:, :, i])
            tor[i, 1:] = om[1:]/(sh.ell[1:]*(sh.ell[1:]+1)) * \
                         gr.radius[i]**2 * rho0[i]

        f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_max))
        tor.byteswap(True)
        tor.tofile(f)
        f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_max))

        tor = np.zeros((gr.n_r_max, gr.lm_max), dtype=np.complex128)
        f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_max))
        tor.byteswap(True)
        tor.tofile(f)
        f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_max))

        # Calculate the pressure
        pre = np.zeros((gr.n_r_max, gr.lm_max), dtype=np.complex128)
        for i in range(gr.n_r_max):
            p = sh.spat_spec(gr.pre[:, :, i])
            pre[i, :] = p[:]

        f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_max))
        pre.byteswap(True)
        pre.tofile(f)
        f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_max))

        pre = np.zeros((gr.n_r_max, gr.lm_max), dtype=np.complex128)
        f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_max))
        pre.byteswap(True)
        pre.tofile(f)
        f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_max))

        # Calculate entropy/temperature
        if l_heat:
            entropy = np.zeros((gr.n_r_max, gr.lm_max), dtype=np.complex128)
            for i in range(gr.n_r_max):
                p = sh.spat_spec(gr.entropy[:, :, i])
                entropy[i, :] = p[:]

            f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_max))
            entropy.byteswap(True)
            entropy.tofile(f)
            f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_max))

            entropy = np.zeros((gr.n_r_max, gr.lm_max), dtype=np.complex128)
            f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_max))
            entropy.byteswap(True)
            entropy.tofile(f)
            f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_max))

        # Calculate chemical composition
        if l_chem:
            xi = np.zeros((gr.n_r_max, gr.lm_max), dtype=np.complex128)
            for i in range(gr.n_r_max):
                p = sh.spat_spec(gr.xi[:, :, i])
                xi[i, :] = p[:]

            f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_max))
            xi.byteswap(True)
            xi.tofile(f)
            f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_max))

            xi = np.zeros((gr.n_r_max, gr.lm_max), dtype=np.complex128)
            f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_max))
            xi.byteswap(True)
            xi.tofile(f)
            f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_max))

        if l_mag:
            # Calculate and store the poloidal potential using Br
            pol = np.zeros((gr.n_r_max, gr.lm_max), dtype=np.complex128)
            for i in range(gr.n_r_max):
                vr = sh.spat_spec(gr.Br[:, :, i])
                pol[i, 1:] = vr[1:]/(sh.ell[1:]*(sh.ell[1:]+1))*gr.radius[i]**2

            f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_max))
            pol.byteswap(True)
            pol.tofile(f)
            f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_max))

            pol = np.zeros((gr.n_r_max, gr.lm_max), dtype=np.complex128)
            f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_max))
            pol.byteswap(True)
            pol.tofile(f)
            f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_max))

            # Calculate the toroidal potential using jr
            tor = np.zeros((gr.n_r_max, gr.lm_max), dtype=np.complex128)

            jr = 1./s3D*(thetaderavg(np.sin(th3D)*gr.Bphi, order=4) -
                         phideravg(gr.Btheta, minc=gr.minc))

            for i in range(gr.n_r_max):
                om = sh.spat_spec(jr[:, :, i])
                tor[i, 1:] = om[1:]/(sh.ell[1:]*(sh.ell[1:]+1))*gr.radius[i]**2

            f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_max))
            tor.byteswap(True)
            tor.tofile(f)
            f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_max))

            tor = np.zeros((gr.n_r_max, gr.lm_max), dtype=np.complex128)
            f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_max))
            tor.byteswap(True)
            tor.tofile(f)
            f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_max))

        if l_mag and l_cond_ic:
            ri = gr.radratio/(1.-gr.radratio)


            # Calculate and store the poloidal potential using Br
            pol = np.zeros((gr.n_r_ic_max, gr.lm_max), dtype=np.complex128)
            for i in range(gr.n_r_ic_max):
                rdep = np.zeros(sh.ell.shape, dtype=np.float64)
                if i == 0:  # ICB radius
                    vr = sh.spat_spec(gr.Br[:, :, -1])
                    rr = gr.radius[-1]
                    rdep = 1.
                else:
                    vr = sh.spat_spec(gr.Br_ic[:, :, i-1])
                    rr = gr.radius_ic[i-1]
                    rdep = (gr.radius_ic[i-1]/ri)**(sh.ell[1:]+1)

                pol[i, 1:] = vr[1:]/(sh.ell[1:]*(sh.ell[1:]+1))*rr**2
                # Technically we should define by rdep but this yields
                # numerical difficulties
                #pol[i, 1:] /= rdep

            f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_ic_max))
            pol.byteswap(True)
            pol.tofile(f)
            f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_ic_max))

            pol = np.zeros((gr.n_r_ic_max, gr.lm_max), dtype=np.complex128)
            f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_ic_max))
            pol.byteswap(True)
            pol.tofile(f)
            f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_ic_max))

            # Calculate the toroidal potential using jr
            tor = np.zeros((gr.n_r_ic_max, gr.lm_max), dtype=np.complex128)

            th3D = np.zeros_like(gr.Br_ic)
            rr3D = np.zeros_like(th3D)
            for i in range(gr.n_theta_max):
                th3D[:, i, :] = gr.colatitude[i]
            for i in range(gr.n_r_ic_max-1):
                rr3D[:, :, i] = gr.radius_ic[i]
            rr3D[:, :, -1] = 1e-4
            s3D = rr3D*np.sin(th3D)
            jr_ic = np.zeros_like(th3D)
            jr_ic = 1./s3D*(thetaderavg(np.sin(th3D)*gr.Bphi_ic, order=4) -
                            phideravg(gr.Btheta_ic, minc=gr.minc))

            for i in range(gr.n_r_ic_max):
                if i == 0:  # ICB radius
                    om = sh.spat_spec(jr[:, :, -1])
                    rr = gr.radius[-1]
                else:
                    om = sh.spat_spec(jr_ic[:, :, i-1])
                    rr = gr.radius_ic[i-1]
                tor[i, 1:] = om[1:]/(sh.ell[1:]*(sh.ell[1:]+1))*rr**2

            f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_ic_max))
            tor.byteswap(True)
            tor.tofile(f)
            f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_ic_max))

            tor = np.zeros((gr.n_r_ic_max, gr.lm_max), dtype=np.complex128)
            f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_ic_max))
            tor.byteswap(True)
            tor.tofile(f)
            f.write(struct.pack('>i', 16*gr.lm_max*gr.n_r_ic_max))

        f.close()


if __name__ == '__main__':

    gr = MagicGraph()
    Graph2Rst(gr, dt=1e-4)
