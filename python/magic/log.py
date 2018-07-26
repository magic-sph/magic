# -*- coding: utf-8 -*-
import re
import os
import numpy as np

class MagicSetup:
    """
    This class allows to read the input namelist or the log file of a current job
    and creates an object that contains all the parameters found in the
    namelist/log file.

    >>> stp = MagicSetup(nml='log.test', quiet=True)
    >>> print(stp.ra) # print the Rayleigh number
    >>> print(stp.n_r_max) # print n_r_max
    """

    def __init__(self, datadir='.', nml='input.nml', quiet=False):
        """
        :param datadir: the working directory
        :type datadir: str
        :param nml: name of the input namelist/ log file
        :type nml: str
        :param quiet: when set to True, makes the output silent (default False)
        :type quiet: bool
        """
        logFile = re.compile(r'log\.(.*)')
        valueInt  = re.compile(r'^([0-9]+)$')
        valueReal = re.compile(r'[+-]?([0-9]+\.[0-9]*|[0-9]*\.[0-9]+)')
        valueNumber = re.compile(r'\b(([\+\-]?[0-9]+)?\.)?[0-9]*([eE][-+]?[0-9]+)?')
        valueFalse = re.compile(r"(\.(true|false|t|f)\.)",re.I)
        valueTrue = re.compile(r"(\.(true|t)\.)",re.I)
        valueNone = re.compile(r"(NONE)",re.I)
        filename = os.path.join(datadir, nml)
        file = open(filename,'r')
        tab = file.readlines()
        tab2 = []
        for i in tab:
            if re.search('=', i) is not None:
                st = i.replace(',', ' ')
                st = st.rstrip('\n')
                if not st.startswith(' !'):
                    tab2.append(st)

        for i in tab2:
            val = i.split('=')
            lhs = val[0].strip()
            lhs = lhs.replace(' ', '_')
            rhs = val[1].strip()
            rhs = rhs.strip('"')
            if valueReal.match(rhs):
                rhs = rhs.replace('D', 'e')
                rhs = rhs.replace('d', 'e')
                try:
                    rhs = float(rhs)
                except ValueError:
                    pass
            elif valueFalse.match(rhs):
                rhs = False
            elif valueTrue.match(rhs):
                rhs = True
            elif valueNone.match(rhs):
                rhs = None
            elif valueInt.match(rhs):
                rhs = int(rhs)
            setattr(self, lhs, rhs)

        self.ra = float(self.ra)
        if not quiet:
            print(self)

        # Overwrite self.tag to be sure that nothing is messed up
        if logFile.match(nml):
            self.tag = logFile.search(nml).groups()[0]

        try:
            # ##############################################
            # Anelastic Benchmark, Jones et al, 2011, Icarus
            # ##############################################
            self.cylinder =  np.arcsin(self.radratio)
            self.WO = ((self.radratio + 1.)
                       /((self.radratio
                          *np.exp(self.strat/self.polind) +1.)))
            self.WI = (self.radratio + 1. - self.WO)/self.radratio
            self.c0 = (2*self.WO - self.radratio -1.)/(1 - self.radratio)
            self.c1 = (1+self.radratio)*(1-self.WO)/(1-self.radratio)**2
            if self.n_lScale==0:
                self.R_out = 1./(1. - self.radratio)
                self.R_inn = self.R_out - 1.
                self.vol = 4./3.*np.pi*(
                    self.R_out**3 - self.R_inn**3)

        except AttributeError as err:
            ## may be not a "standard" anelastic run
            ## check default values or look at the log.TAG file
            if not quiet:
                print(err)

    #def __repr__(self):
        #st = 'Welcome in the run %s\n' % self.tag
        #st += ' ---- Params ---- \n'
        #st += 'Rayleigh = %.2e\n' % self.ra
        #st += 'Ekman = %.2e\n' % self.ek
        #st += 'Prandtl = %.2e\n' % self.pr
        #st += 'Magnetic Prandtl = %.2e' % self.prmag
        #return st
