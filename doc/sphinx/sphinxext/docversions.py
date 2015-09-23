# -*- coding: utf8 -*-
"""Sphinx directive to list other available versions of the current documentation"""
# Copyright or © or Copr. Actimar/IFREMER (2014-2015)
#
# This software is a computer program whose purpose is to provide
# utilities for handling oceanographic and atmospheric data,
# with the ultimate goal of validating the MARS model from IFREMER.
#
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
#
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
#
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
#
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.

import os, inspect, sys, re
from collections import OrderedDict
from glob import glob

from sphinx.directives import Directive
from sphinx.util.console import bold
from docutils.parsers.rst.directives import unchanged,single_char_or_unicode,positive_int
from docutils import nodes
from docutils.statemachine import string2lines
class VersionFinder(object):

    def __init__(self, config):

        # Store config
        self.subpath_doc = config.docversions_subpath_doc.strip(os.path.sep)
        self.subpath_buildhtml = os.path.join(self.subpath_doc,
            config.docversions_subpath_buildhtml)
        self.subpath_version_file = config.docversions_subpath_version_file
        self.subpath_version_var = config.docversions_subpath_version_var
        self.subpath_versions = config.docversions_subpath_versions
        self.index_html = config.docversions_index_html
        self.template = config.docversions_template

        # Where are we?
        self.current_src_dir = os.path.abspath(
            os.path.join(*['..']*(self.subpath_doc.count(os.path.sep)+1)))
        self.current_html_dir = os.path.join(self.current_src_dir, self.subpath_buildhtml)
        if os.path.basename(self.current_src_dir) =='trunk' :
            self.trunk_path = self.current_src_dir
            self.tags_path = os.path.abspath(os.path.join(self.trunk_path, '../tags'))
            self.branches_path = os.path.abspath(
                os.path.join(self.trunk_path, '../branches'))
        else:
            self.trunk_path = os.path.abspath(
                os.path.join(self.current_src_dir, '../../trunk'))
            upper_path = os.path.abspath(
                os.path.join(self.current_src_dir, '..'))
            if upper_path=='tags':
                self.tags_path = upper_path
                self.branches_path = os.path.abspath(
                    os.path.join(self.current_src_dir, '../../branches'))
            else:
                self.branches_path = upper_path
                self.tags_path = os.path.abspath(
                    os.path.join(self.current_src_dir, '../../tags'))


    def find_version_label(self, src_dir):
        """Find the package version knowing the src dir"""
        version_file = os.path.join(src_dir, self.subpath_version_file)
        if not os.path.exists(version_file): return
        try:
            f = open(version_file)
            for line in f:
                if line.startswith(self.subpath_version_var):
                    exec(line[:-1])
                    f.close()
                    return eval(self.subpath_version_var)
            f.close()
        except Exception as e:
            print 'Cannot retreive version label for dir: '+src_dir
            print e.message


    def get_tag_specs(self, src_dir):
        """Get the version name and directories knowing the src dir"""
        if not os.path.isdir(src_dir): return
        if src_dir==self.current_src_dir: return

        # Version label
        version = self.find_version_label(src_dir)
        if not version: return

        # Paths
        html_dir = os.path.join(src_dir, self.subpath_buildhtml)
        if os.path.isdir(html_dir):
            return (version, html_dir)

    def get_branch_specs(self, src_dir):
        """Get the branch name, version name and directories knowing the src dir"""
        specs = self.get_tag_specs(src_dir)
        if specs is None: return
        return os.path.basename(src_dir), specs[1]




class DocversionsDirective(Directive):

    has_content = False
    required_arguments = 0
    optional_arguments = 0
    option_spec = {}

    def run(self):
        """Generate the list and make symlinks"""
        # Init
        env = self.state.document.settings.env
        versions = env.docversions_list
        current_html_dir = env.docversions_current_html_dir
        index = env.docversions_index_html
        verdir = env.docversions_subpath_versions
        all_dir = os.path.join(current_html_dir, verdir)
        if not os.path.exists(all_dir):
            os.makedirs(all_dir)

        # Links + rst
        rst_links = []
        for name, html_dir in versions:
            symlink = os.path.join(all_dir, name)
            if os.path.islink(symlink):
                os.remove(symlink)
            os.symlink(html_dir, os.path.join(all_dir, name))
            rst_links.append('`%(name)s <%(verdir)s/%(name)s/%(index)s>`_'%locals())
        versions = ', '.join(rst_links)
        raw_text = env.docversions_template%locals()

        # Insert
        source = self.state_machine.input_lines.source(self.lineno - self.state_machine.input_offset - 1)
        include_lines = string2lines(raw_text, convert_whitespace=1)
        self.state_machine.insert_input(include_lines,source)

        return []

def store_versions(app):
    """Find versions and store them in app.builder.env"""

    app.info(bold('searching for other html docs...'), nonl=True)

    # Version finder
    VF = VersionFinder(app.config)

    # Find directory and name of versions
    versions = []
    # - trunk
    specs = VF.get_tag_specs(VF.trunk_path)
    if specs:
        versions.append(('trunk',)+specs[1:])
    # - tags
    if os.path.isdir(VF.tags_path):
        for src_dir in sorted(glob(VF.tags_path+'/*')):
            specs = VF.get_tag_specs(src_dir)
            if specs:
                versions.append(specs)
    # - branches
    if os.path.isdir(VF.branches_path):
        for src_dir in sorted(glob(VF.branches_path+'/*')):
            specs = VF.get_branch_specs(src_dir)
            if specs:
                versions.append(specs)


    # Store them
    app.builder.env.docversions_list = versions
    app.builder.env.docversions_current_html_dir = VF.current_html_dir
    app.builder.env.docversions_subpath_versions = VF.subpath_versions
    app.builder.env.docversions_index_html = VF.index_html
    app.builder.env.docversions_template = VF.template

    app.info('done')

def setup(app):
    app.add_config_value('docversions_subpath_doc', 'doc', 'html')
    app.add_config_value('docversions_subpath_buildhtml', 'build/html', 'html')
    app.add_config_value('docversions_index_html', 'index.html', 'html')
    app.add_config_value('docversions_subpath_version_file', 'lib/python/vacumm/__init__.py', 'html')
    app.add_config_value('docversions_subpath_version_var', '__version__', 'html')
    app.add_config_value('docversions_subpath_versions', 'versions', 'html')
    app.add_config_value('docversions_template', 'Other versions: %(versions)s', 'html')
    app.connect('builder-inited', store_versions)
    app.add_directive('docversions', DocversionsDirective)

