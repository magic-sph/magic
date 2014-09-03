# -*- coding: utf-8 -*-
import re
import os
import string

__author__  = "$Author: gastine $"
__date__   = "$Date: 2013-12-02 16:36:18 +0100 (Mon, 02 Dec 2013) $"
__version__ = "$Revision: 553 $"


class MagicSetup:

    def __init__(self, datadir='.', nml='input.nml', quiet=False):
        logFile = re.compile(r'log\.(.*)')
        valueInt  = re.compile(r'^[0-9]$')
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
                st = string.replace(i, ',', ' ')
                st = st.rstrip('\n')
                if not st.startswith(' !'):
                    tab2.append(st)

        for i in tab2:
            val = string.split(i, '=')
            lhs = val[0].strip()
            lhs = lhs.replace(' ', '_')
            rhs = val[1].strip()
            rhs = rhs.strip('"')
            if valueReal.match(rhs):
                rhs = string.replace(rhs, 'D', 'e')
                rhs = string.replace(rhs, 'd', 'e')
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
            print self

        # Overwrite self.tag to be sure that nothing is messed up
        if logFile.match(nml):
            self.tag = logFile.search(nml).groups()[0]


    #def __repr__(self):
        #st = 'Welcome in the run %s\n' % self.tag
        #st += ' ---- Params ---- \n'
        #st += 'Rayleigh = %.2e\n' % self.ra
        #st += 'Ekman = %.2e\n' % self.ek
        #st += 'Prandtl = %.2e\n' % self.pr
        #st += 'Magnetic Prandtl = %.2e' % self.prmag
        #return st
