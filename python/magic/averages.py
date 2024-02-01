#,-*- coding: utf-8 -*-
import os
import json
from types import SimpleNamespace
from collections import OrderedDict
import numpy as np
from .libmagic import scanDir, avgField
from .series import MagicTs
from .log import MagicSetup
from .spectrum import MagicSpectrum
from .radial import MagicRadial
from .libmelt import MagicMelt
from .movie import Movie

if 'MAGIC_HOME' in os.environ:
    default_model = os.path.join(os.environ['MAGIC_HOME'], 'python/magic/model.json')
else:
    default_model = os.path.join(os.environ['HOME'], 'magic/python/magic/model.json')

INDENT = 3
SPACE = " "
NEWLINE = "\n"
def to_json(o, level=0):
    """
    This is a manual JSON serializer function that makes the outputs look a little
    bit better than default.
    """
    ret = ""
    if isinstance(o, dict):
        ret += "{" + NEWLINE
        comma = ""
        for k, v in o.items():
            ret += comma
            comma = ",\n"
            ret += SPACE * INDENT * (level + 1)
            ret += '"' + str(k) + '":' + SPACE
            ret += to_json(v, level + 1)
  
        ret += NEWLINE + SPACE * INDENT * level + "}"
    elif isinstance(o, str):
        ret += '"' + o + '"'
    elif isinstance(o, list):
        ret += "[" + ",".join([to_json(e, level + 1) for e in o]) + "]"
    # Tuples are interpreted as lists
    elif isinstance(o, tuple):
        ret += "[" + ",".join(to_json(e, level + 1) for e in o) + "]"
    elif isinstance(o, bool):
        ret += "true" if o else "false"
    elif isinstance(o, int):
        ret += str(o)
    elif isinstance(o, float):
        if abs(o) > 1e2 or abs(o) < 1e-2:
            ret += '{:.8e}'.format(o)
        else:
            ret += '{:.8g}'.format(o)
    elif isinstance(o, np.ndarray) and np.issubdtype(o.dtype, np.integer):
        ret += "[" + ','.join(map(str, o.flatten().tolist())) + "]"
    elif isinstance(o, np.ndarray) and np.issubdtype(o.dtype, np.bool_):
        ret += "[" + ','.join(map(lambda x: '"{}"'.format(x), o.flatten().tolist())) + "]"
    elif isinstance(o, np.ndarray) and np.issubdtype(o.dtype, np.str):
        ret += "[" + ','.join(map(lambda x: '"{}"'.format(x), o.flatten().tolist())) + "]"
    elif isinstance(o, np.ndarray) and np.issubdtype(o.dtype, np.inexact):
        ret += "[" + ','.join(map(lambda x: '{:.8e}'.format(x), o.flatten().tolist())) + "]"
    elif o is None:
        ret += 'null'
    else:
        raise TypeError("Unknown type '{}' for json serialization".format(str(type(o))))

    return ret


class AvgField:
    """
    This class computes the time-average properties from time series, spectra
    and radial profiles. It will store the input starting time in a small file
    named ``tInitAvg``, such that the next time you use it you don't need to
    provide ``tstart`` again. By default, the outputs are stored in a 
    fully documented JSON file named avg.json: this is split into several
    categories, namely numerical parameters, physical parameters, time averaged
    scalar quantities, time averaged spectra and time-averaged radial profiles.
    The quantities stored in the JSON file are entirely controlled by an input
    model file wich enlists the different quantities of interest. 
    A default example file named ``model.json`` is provided in
    ``$MAGIC_HOME/python/magic``, and an example of how to build a dedicated one
    is provided below.

    >>> # Average from t=2.11
    >>> a = AvgField(tstart=2.11)
    >>> # Average only the files that match the pattern N0m2[a-c]
    >>> a = AvgField(tstart=2.11, tag='N0m2[a-c]')
    >>> print(a) # print the formatted output
    >>> # Custom JSON model to select averages
    >>> json_model = { 'phys_params': ['ek'],
                       'time_series': { 'heat': ['topnuss', 'botnuss'],
                                        'e_kin': ['ekin_pol', 'ekin_tor'],
                                        'par': ['rm'] },
                       'spectra': {},
                       'radial_profiles': {'powerR': ['viscDiss', 'buoPower']}
                     }
    >>> # Compute the selected averages in the dirctory mydir                 
    >>> a = AvgField(datadir='mydir', model=json_model)
    """

    def __init__(self, tag=None, tstart=None, model=default_model,
                 datadir='.', std=False, write=True):
        """
        :param tag: if you specify an input tag (generic regExp pattern),
                    the averaging process will only happen on the time series
                    that match this input pattern
        :type tag: str
        :param tstart: the starting time for averaging
        :type tstart: float
        :param datadir: working directory
        :type datadir: str
        :param std: compute the standard deviation when set to True
        :type std: bool
        :param write: write the outputs in a JSON file
        :type write: bool
        :param model: this is the path of a JSON file which defines which fields will be
                      handled in the time-averaging process. This can be any
                      python attributes wich is defined in MagicTs, MagicSpectrum
                      or MagicRadial.
        :type model: str
        """

        if not os.path.exists(datadir):
            print('Directory "{}" has not been found'.format(datadir))
            return

        tInitFile = os.path.join(datadir, 'tInitAvg')
        if os.path.exists(tInitFile) and tstart is None:
            with open(tInitFile, 'r') as f:
                st = f.readline().strip('\n')
                tstart = float(st)
        elif tstart is not None:
            with open(tInitFile, 'w') as f:
                f.write('{}'.format(tstart))

        if type(model) == str:
            with open(model, 'r') as f:
                params = json.load(f)
        else: # This is directly a json dict
            params = model

        pattern = os.path.join(datadir, 'log.*')
        logs = scanDir(pattern)

        self.lut = OrderedDict()
        if 'phys_params' in params:
            self.lut['phys_params'] = {}
        if 'num_params' in params:
            self.lut['num_params'] = {}

        # First grab the requested control parameters
        if len(logs) > 0:
            stp = MagicSetup(nml=logs[-1], quiet=True)
            if 'phys_params' in self.lut:
                for p in params['phys_params']:
                    if hasattr(stp, p):
                        self.lut['phys_params'][p] = getattr(stp, p)
                        setattr(self, p, getattr(stp, p))
            if 'num_params' in self.lut:
                for p in params['num_params']:
                    if hasattr(stp, p):
                        self.lut['num_params'][p] = getattr(stp, p)
                        setattr(self, p, getattr(stp, p))

        # Handle time series
        self.lut['time_series'] = {}
        for k, key in enumerate(params['time_series'].keys()):
            ts = MagicTs(field=key, datadir=datadir, all=True, tag=tag, iplot=False)
            # Number of columns has been possibly modified for those files
            if key == 'par' or key == 'geos' or key == 'heat' or key == 'dtVrms' \
               or key == 'phase':
                fix = True
            else:
                fix = False
            if hasattr(ts, 'time'): # Manage to read file
                mask = np.where(abs(ts.time-tstart) == min(abs(ts.time-tstart)), 1, 0)
                ind = np.nonzero(mask)[0][0]

                if k == 0:
                    self.lut['time_series']['tavg'] = ts.time[-1]-ts.time[ind]
                    self.lut['time_series']['total_time'] = ts.time[-1]-ts.time[0]
                    self.tavg = ts.time[-1]-ts.time[ind]
                    self.trun = ts.time[-1]-ts.time[0]

                for field in params['time_series'][key]:
                    if hasattr(ts, field):
                        if std and field != 'dt':
                            xmean, xstd = avgField(ts.time[ind:],
                                                   ts.__dict__[field][ind:],
                                                   std=True,
                                                   fix_missing_series=fix)
                            self.lut['time_series'][field+'_av'] = xmean
                            self.lut['time_series'][field+'_sd'] = xstd
                            setattr(self, field+'_av', xmean)
                            setattr(self, field+'_sd', xstd)
                        else:
                            xmean = avgField(ts.time[ind:],
                                             ts.__dict__[field][ind:],
                                             fix_missing_series=fix)
                            self.lut['time_series'][field+'_av'] = xmean
                            setattr(self, field+'_av', xmean)
            else:  # If parameters is absent then put it to -1
                for field in params['time_series'][key]:
                    self.lut['time_series'][field+'_av'] = -1
                    setattr(self, field+'_av', -1)
                    if std:
                        self.lut['time_series'][field+'_sd'] = -1
                        setattr(self, field+'_sd', -1)


        # Get tags involved in averaging for spectra and radial profiles
        tags = self.get_tags(tstart)

        # Handle spectra
        self.lut['spectra'] = {}
        # Determine whether file exists
        if len(tags) > 0:
            file_exists = True
            # If only one tag is retained but averaged file does not exist yet
            if len(tags) == 1:
                file = os.path.join(datadir, 'kin_spec_ave.' + tags[-1])
                if not os.path.exists(file):
                    file_exists = False
        else:
            file_exists = False

        if file_exists:
            for key in params['spectra'].keys():
                sp = MagicSpectrum(field=key, datadir=datadir, iplot=False, tags=tags,
                                   quiet=True)
                if hasattr(sp, 'index'): # Manage to read file
                    self.lut['spectra']['index'] = sp.index
                    for field in params['spectra'][key]:
                        if hasattr(sp, field):
                            self.lut['spectra'][field+'_spec_av'] = sp.__dict__[field]
                            if std and hasattr(sp, field + '_SD'):
                                self.lut['spectra'][field+'_spec_sd'] = \
                                    sp.__dict__[field + '_SD']
        else:  # Set parameters to -1
            for key in params['spectra'].keys():
                self.lut['spectra']['index'] = -1 * np.ones(32)
                for field in params['spectra'][key]:
                    self.lut['spectra'][field+'_spec_av'] = -1 * np.ones(32)
                    if std:
                        self.lut['spectra'][field+'_spec_sd'] = -1 * np.ones(32)

        # Handle radial profiles
        self.lut['radial_profiles'] = {}
        if len(tags) > 0:
            for key in params['radial_profiles'].keys():
                rr = MagicRadial(field=key, datadir=datadir, iplot=False, tags=tags,
                                 quiet=True)
                if hasattr(rr, 'radius'): # Manage to read file
                    self.lut['radial_profiles']['radius'] = rr.radius
                    for field in params['radial_profiles'][key]:
                        if hasattr(rr, field):
                            self.lut['radial_profiles'][field+'_av'] = rr.__dict__[field]
                            setattr(self, field+'R_av', rr.__dict__[field])
                            if std and hasattr(rr, field + '_SD'):
                                self.lut['radial_profiles'][field+'_sd'] = \
                                    rr.__dict__[field + '_SD']
                                setattr(self, field+'R_sd', rr.__dict__[field+'_SD'])
        else:  # Set parameters to -1
            for key in params['radial_profiles'].keys():
                self.lut['radial_profiles']['radius'] = -1 * np.ones(33)
                for field in params['radial_profiles'][key]:
                    self.lut['radial_profiles'][field+'_av'] = -1 * np.ones(33)
                    #setattr(self, field+'R_av', rr.__dict__[field])
                    if std:
                        self.lut['radial_profiles'][field+'_sd'] = -1 * np.ones(33)
                        #setattr(self, field+'R_sd', rr.__dict__[field+'_SD'])

        # Handle theta profiles
        self.lut['theta_profiles'] = {}
        if 'theta_profiles' in params.keys():
            for key in params['theta_profiles']:

                if key == 'rmelt_theta':
                    files = scanDir(os.path.join(datadir, 'rmelt.*'))
                    if len(files) != 0:
                        rm = MagicMelt(all=True, datadir=datadir, iplot=False)
                        if hasattr(rm, 'colatitude'):
                            self.lut['theta_profiles']['colatitude'] = rm.colatitude
                            mask = np.where(abs(rm.time-tstart) == \
                                            min(abs(rm.time-tstart)), 1, 0)
                            ind = np.nonzero(mask)[0][0]

                            if std:
                                xmean, xstd = avgField(rm.time[ind:],
                                                       rm.rmelt[ind:, :],
                                                       std=True)
                                self.lut['theta_profiles']['rmelt_theta_av'] = xmean
                                self.lut['theta_profiles']['rmelt_theta_sd'] = xstd
                                setattr(self, 'rmelt_theta_av', xmean)
                                setattr(self, 'rmelt_theta_sd', xstd)
                            else:
                                xmean, xstd = avgField(rm.time[ind:],
                                                       rm.rmelt[ind:, :])
                                self.lut['theta_profiles']['rmelt_theta_av'] = xmean
                                setattr(self, 'rmelt_theta_av', xmean)
                elif key == 'heatf_theta':
                    files = scanDir(os.path.join(datadir, 'AHF_mov.*'))
                    if len(files) != 0:
                        for k, file in enumerate(files):
                            if k == 0:
                                m = Movie(file=file, iplot=False, datadir=datadir)
                            else:
                                m += Movie(file=file, iplot=False, datadir=datadir)
                        if hasattr(m, 'theta'):
                            data = m.data[0, :, :, 0]  # dT/dr at outer boundary
                            self.lut['theta_profiles']['colatitude'] = m.theta
                            mask = np.where(abs(m.time-tstart) == \
                                            min(abs(m.time-tstart)), 1, 0)
                            ind = np.nonzero(mask)[0][0]

                            if std:
                                xmean, xstd = avgField(m.time[ind:],
                                                       data[ind:, :],
                                                       std=True)
                                self.lut['theta_profiles']['heatf_theta_av'] = xmean
                                self.lut['theta_profiles']['heatf_theta_sd'] = xstd
                                setattr(self, 'heatf_theta_av', xmean)
                                setattr(self, 'heatf_theta_sd', xstd)
                            else:
                                xmean, xstd = avgField(m.time[ind:],
                                                       data[ind:, :])
                                self.lut['theta_profiles']['heatf_theta_av'] = xmean
                                setattr(self, 'heatf_theta_av', xmean)

        # Write a json file
        if write:
            self.write_json(datadir)

    def write_header(self):
        """
        Write header in case an ascii output is requested
        """
        st = '#'

        if 'phys_params' in self.lut:
            for key in self.lut['phys_params']:
                st += key.rjust(max(10, len(key)+1))
        for key in self.lut['time_series']:
            st += key.rjust(max(16, len(key)+1))
        if 'num_params' in self.lut:
            for key in self.lut['num_params']:
                par = self.lut['num_params'][key]
                if type(par) != str and type(par) != bool:
                    st += key.rjust(len(key)+1)

        return st

    def __str__(self):
        """
        Formatted output
        """
        st = ' '

        if 'phys_params' in self.lut:
            for par in self.lut['phys_params'].values():
                st += '{:10.3e}'.format(par)
        for par in self.lut['time_series'].values():
            st += '{:16.8e}'.format(par)
        if 'num_params' in self.lut:
            for par in self.lut['num_params'].values():
                if type(par) != str and type(par) != bool:
                    st += ' {:g}'.format(par)

        return st

    def get_tags(self, tstart):
        """
        This routine returns a list of tags which have been generated after tstart

        :param tstart: starting averaging time
        :type tstart: float
        :returns: a list of tags
        :rtype: list
        """
        logFiles = scanDir('log.*')
        tags = []
        for lg in logFiles:
            nml = MagicSetup(nml=lg, quiet=True)
            if nml.start_time > tstart:
                tags.append(nml.tag)

        return tags

    def write_json(self, datadir='.'):
        """
        This function writes the averages as a simple JSON file stored in the
        directory 'avg.json'

        :param datadir: working directory
        :type datadir: str
        """
        with open(os.path.join(datadir, 'avg.json'), 'w') as f:
            st = to_json(self.lut)
            f.write(st)


class AvgStack:
    """
    This class is made to go through a list of directories and gather all the local
    avg files to compute a global output which summarises all the outputs in one
    single file

    >>> # Simply go through the directories listed in "runs.txt" and produce a 
    >>> # local file named "my_results.json"
    >>> st = AvgStack(dirList="runs.txt", filename="my_results.json")
    >>> # Now also recompute each individual avg.json file in each directory
    >>> st = AvgStack(dirList="runs.txt", filename="my_results.json",
                      recompute=True, module="path_to_model.json", std=True)
    """

    def __init__(self, dirList='runs.txt', filename='avg.json', datadir='.',
                 recompute=False, model=default_model, std=False, readonly=False):
        """
        :param dirList: the filename of the list of directories, lines starting with
                        a hash are omitted by the reader
        :type dirList: str
        :param filename: 
        :param recompute: recompute each individual average file in each single 
                          directory when set to True
        :type recompute: bool
        :param datadir: working directory
        :type datadir: str
        :param std: compute the standard deviation when set to True
        :type std: bool
        :param model: this is the path of a JSON file which defines which fields will be
                      handled in the time-averaging process. This can be any
                      python attributes wich is defined in MagicTs, MagicSpectrum
                      or MagicRadial.
        :type model: str
        :param readonly: when set to True, simply read the summary json file
        :param readonly: bool
        """

        if readonly: # Only read the json file
            with open(os.path.join(datadir, filename), 'r') as f:
                #self.lut = json.load(f, object_hook=lambda d: SimpleNamespace(**d))
                self.lut = json.load(f)
        else:
            with open(os.path.join(datadir, dirList), 'r') as f:
                dirs = f.readlines()

            startdir = os.getcwd()
            for k, dir in enumerate(dirs):
                if not dir.startswith('#'):
                    print("In dir {}".format(dir.rstrip('\n')))
                    os.chdir(os.path.join(datadir, dir.rstrip('\n')))
                    if recompute:
                        a = AvgField(model=model, std=std)
                    if not hasattr(self, 'lut'):
                        self.lut = self.load()
                        keys = ['spectra', 'radial_profiles', 'theta_profiles']
                        for key in keys:
                            if key in self.lut.keys():
                                for key1 in self.lut[key].keys():
                                    self.lut[key][key1] = [np.atleast_1d(self.lut[key][key1])]
                    else:
                        lut = self.load()
                        self.merge(lut)

                    os.chdir(startdir)

            with open(os.path.join(datadir, filename), 'w') as f:
                st = to_json(self.lut)
                f.write(st)

        self.simple_namespace()

    def __add__(self, new):
        """
        This routine is used to add two AvgStack objects.

        :param new: the lookup table which needs to be added
        :type new: AvgStack
        """
        keys = ['phys_params', 'num_params', 'time_series']
        for key in keys:
            if key in new.lut.keys() and key in self.lut.keys():
                for key1 in self.lut[key].keys():
                    if key1 in new.lut[key].keys():
                        arr = np.atleast_1d(self.lut[key][key1])
                        arr1 = np.atleast_1d(new.lut[key][key1])
                        self.lut[key][key1] = np.concatenate((arr, arr1))

        keys = ['spectra', 'radial_profiles', 'theta_profiles']
        for key in keys:
            if key in new.lut.keys() and key in self.lut.keys():
                for key1 in self.lut[key].keys():
                    if key1 in new.lut[key].keys():
                        for lst in new.lut[key][key1]:
                            arr1 = np.atleast_1d(lst)
                            self.lut[key][key1].append(arr1)

        self.simple_namespace()

        return self

    def simple_namespace(self):
        """
        This routine creates a simpler namespace from the lookup table
        """
        keys = ['phys_params', 'num_params', 'time_series']
        for key in keys:
            if key in self.lut.keys():
                for key1 in self.lut[key].keys():
                    setattr(self, key1, np.atleast_1d(self.lut[key][key1]))

        keys = ['spectra', 'radial_profiles', 'theta_profiles']
        for key in keys:
            if key in self.lut.keys():
                for key1 in self.lut[key].keys():
                    g = []
                    for ff in self.lut[key][key1]:
                        g.append(np.atleast_1d(ff))
                    setattr(self, key1, g)

    def merge(self, newlut):
        """
        This routine is used to merge two lookup tables.

        :param newlut: the lookup table which needs to be added
        :type newlut: dict
        """
        keys = ['phys_params', 'num_params', 'time_series']
        for key in keys:
            if key in newlut.keys() and key in self.lut.keys():
                for key1 in self.lut[key].keys():
                    if key1 in newlut[key].keys():
                        arr = np.atleast_1d(self.lut[key][key1])
                        arr1 = np.atleast_1d(newlut[key][key1])
                        self.lut[key][key1] = np.concatenate((arr, arr1))

        keys = ['spectra', 'radial_profiles', 'theta_profiles']
        for key in keys:
            if key in newlut.keys() and key in self.lut.keys():
                for key1 in self.lut[key].keys():
                    if key1 in newlut[key].keys():
                        arr1 = np.atleast_1d(newlut[key][key1])
                        self.lut[key][key1].append(arr1)
                    else: #  Fill with -1 dummy arrays to make sure everything has the right size
                        arr1 = -1 * np.ones(32)
                        self.lut[key][key1].append(arr1)

    def load(self):
        """
        This routine reads the avg.json file stored in one local directory which
        contains a numerical simulation. It returns the corresponding lookup table.

        :returns: a lookup table
        :rtype: dict
        """
        with open('avg.json', 'r') as f:
            lut = json.load(f)

        return lut


if __name__ == '__main__':
    a = Avg(std=True)
    st = a.write_header()
    print(st)
    print(a)
