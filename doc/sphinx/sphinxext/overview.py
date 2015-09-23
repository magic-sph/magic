# -*- coding: utf8 -*-
"""Sphinx directive to add an overview of a python module or class"""
# Copyright or © or Copr. Actimar/IFREMER (2010-2015)
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
#
from sphinx.directives import Directive
from docutils.parsers.rst.directives import unchanged,single_char_or_unicode,positive_int
from docutils import nodes
from docutils.statemachine import string2lines
import inspect, sys, re




def setup(app):
    app.add_config_value('overview_underline', '-',  False)
    app.add_config_value('overview_title_overview', 'Overview', False)
    app.add_config_value('overview_title_content', 'Content', False)
    app.add_config_value('overview_columns', 3, False)

    app.add_directive('overview', OverViewDirective)


class overview(nodes.General, nodes.Element):
    pass

def overview_strings(arg):
    return re.split('[\s,]+', arg.strip())

class OverViewDirective(Directive):
    has_content = True
    option_spec = {}
    option_spec['underline'] = single_char_or_unicode
    option_spec['title-overview'] = unchanged
    option_spec['title-content'] = unchanged
    option_spec['extra-attributes'] = overview_strings
    option_spec['extra-functions'] = overview_strings
    option_spec['extra-class-attributes'] = overview_strings
    option_spec['extra-classes'] = overview_strings
    option_spec['extra-methods'] = overview_strings
    option_spec['columns'] = positive_int
    option_spec['inherited-members'] = unchanged
    required_arguments = 1
    optional_arguments = 0

    def run(self):

        # Get object
        objname = self.arguments[0]
        try:
            __import__(objname)
            object = sys.modules[objname]
        except:
            self.warning('Cannot import object %s for overview'%objname)
            return []

        # Options
        config = self.state.document.settings.env.config
        # - titles
        titles = {}
        for title_name in 'overview', 'content':
            if self.options.has_key('title_'+title_name) and self.options['title_'+title_name] is not None:
                title = self.options['title_'+title_name]
            else:
                title = getattr(config, 'overview_title_'+title_name)
            if not isinstance(title, (str, unicode)) or not title:
                title = False
            titles[title_name] = title
        # - underline
        if self.options.has_key('underline') and self.options['underline'] is not None:
            underline = self.options['underline']
        else:
            underline = config.overview_underline
        underline = str(underline)[0]
        # - extra
        extra={}
        for etype in 'attributes', 'functions', 'classes', 'methods', 'class_attributes':
            etype = 'extra_'+etype
            if self.options.has_key(etype) and self.options[etype] is not None:
                extra[etype] = self.options[etype]
        # - columns
        columns = self.options.get('columns', config.overview_columns)
        # - inheritance
        extra['inherited'] = 'inherited-members' in self.options

        # Format
        raw_text = OverView(object, **extra).format(indent=0,
            title_overview=titles['overview'], title_content=titles['content'],
            underline=underline, columns=columns)
        source = self.state_machine.input_lines.source(self.lineno - self.state_machine.input_offset - 1)
        include_lines = string2lines(raw_text, convert_whitespace=1)
        self.state_machine.insert_input(include_lines,source)


        return []


class OverView(object):
    """Python object rst overview generator

    :Usage:

    >>> import mymodule
    >>> rst_text = OverView(mymodule).format()
    """
    columns = 3
    def __init__(self, object, extra_attributes=[], extra_functions=[], extra_classes=[],
                 extra_methods=[], extra_class_attributes=[], inherited=True):
        self.inherited = inherited

        # Check must be a module or a class
        if not inspect.ismodule(object) and not inspect.isclass(object): raise
        if inspect.ismodule(object):
            self.module = object
        else:
            self.module = inspect.getmodule(object)
        self.modname = self.module.__name__

        # Get base lists
        self.attributes = self.get_members(object)
        self.functions = self.get_members(object, 'function')
        self.classes = self.get_members(object, 'class')

        # Sub content
        self.class_contents = {}
        for clsname, cls in self.classes:
            self.class_contents[clsname] = dict(
                methods=self.get_members(cls, 'method'),
                #classmethods = self.get_members(object, 'classmethod'),
                #staticmethods = self.get_members(object, 'staticmethod'),
                attributes=self.get_members(cls))


        # Check extra args
        extra_names = 'attributes', 'functions', 'classes', 'methods', 'class_attributes'
        for etype in extra_names:
            # Get
            etype = 'extra_'+etype
            extras = eval(etype)
            # Store
            setattr(self, etype, extras)
            # Check prefix
            for i, extra in enumerate(extras):
                if not extra.startswith(self.modname+'.'):
                   extras[i] = self.modname+'.'+extra
        if len(extra_methods+extra_attributes): # Check missing classes
            for objname in extra_methods+extra_attributes:
                clsname = objname.split('.')[1]
                if clsname not in extra_classes:
                    extra_classes.append(clsname)


    def get_extra_class_attributes(self, clsname):
        """Get the list of extra attributes that belongs to a class"""
        modname = self.modname
        return [attr for attr in self.extra_class_attributes if attr.startswith('%(modname)s.%(clsname)s.'%locals())]

    def get_extra_methods(self, clsname):
        """Get the list of extra methods that belongs to a class"""
        modname = self.modname
        return [meth for meth in self.extra_methods if attr.startswith('%(modname)s.%(clsname)s.'%locals())]

    def get_members(self, object, predicate=None):
        """Get the list of object members of a given type"""

        # Get base list
        predicate_spec = predicate
        if predicate is not None:
            if 'is'+predicate in dir(inspect):
                ismatched = predicate = getattr(inspect, 'is'+predicate)
            else:
                ismatched = predicate =  lambda o: isinstance(o, eval(predicate_spec))
        else:
            ismatched = lambda o: True
        if hasattr(object, '__all__'): # Fixed list

            members = [(mname, getattr(object, mname)) for mname in object.__all__
                if (hasattr(object, mname) and ismatched(getattr(object, mname)))]

        else: # Auto list

            # All members
            members = [(mname, member) for mname, member in inspect.getmembers(object, predicate) if not mname.startswith('_')]

             # Inheritance
            if self.inherited is False and hasattr(object, '__dict__'):
                members = [(mname, member) for mname, member in members if mname in object.__dict__.keys()]

            # Filter out non local members
            if not self.inherited or predicate_spec not in ['method', None, 'classmethod', 'staticmethod']:
                members = [(mname, member) for mname, member in members
                    if inspect.getmodule(member) is None or inspect.getmodule(member) is self.module]


        # Attributes only
        if predicate is None:
            members = [(mname, member) for mname, member in members if not inspect.ismethod(member) and
                not inspect.isclass(member) and not inspect.isfunction(member)]
        return members


    @classmethod
    def indent(cls, indent, *text, **kwargs):
        xindent = kwargs.get('xindent', '')
        return '\n'.join([(indent*'\t'+xindent+line) for line in text])

    def format_ref(self, objname, object, clsname=None):
        """Format a reference link to an object"""
        # Declaration type
        if inspect.isfunction(object):
            dectype = 'func'
        elif inspect.isclass(object):
            dectype = 'class'
        elif inspect.ismethod(object) and not isinstance(object, property):
            dectype = 'meth'
        else:
            dectype = 'attr'

        # Class content
        if clsname is None and hasattr(object, 'im_class'):
            clsname = object.im_class.__name__
        if clsname is not None: #FIXME: properties
            objname = '%s.%s'%(clsname, objname)

        # Format
        modname = self.modname
        rst = ":%(dectype)s:`~%(modname)s.%(objname)s`"%locals()
        return rst

    def format_list(self, args, indent=0, columns=None, xindent=''):
        if len(args) == 0: return ''
        if columns is None: columns = self.columns
        if columns == 0:# or len(args) <= columns:
            return ', '.join(args)
        rst = '\n'
        rst += self.indent(indent+1, '.. hlist::\n', xindent=xindent)
        rst += self.indent(indent+2, ':columns: %i\n\n'%columns, xindent=xindent)
        for arg in args:
            rst += self.indent(indent+2,'- %s\n'%arg, xindent=xindent)
        rst += '\n'
        return rst



    def format_title(self, title, underline, indent=0):
        """Format the title of the overview or content paragraphs"""
        rst = ''
        if title:
            rst += self.indent(indent, title, len(title)*underline)
            rst +='\n\n'
        return rst

    def format_attributes(self, indent=0, columns=None):
        """Format module level attributes"""
        rst = ''
        if len(self.attributes)+len(self.extra_attributes):
            rst += self.indent(indent, ':Attributes: ')
            attrs = []
            if self.attributes:
                attrs += [self.format_ref(attname, att) for attname, att in self.attributes]
                #rst += self.format_list(attrs, indent=indent, columns=columns)
                #rst +=', '.join([self.format_ref(attname, att) for attname, att in self.attributes])
            if self.extra_attributes:
                attrs += [':attr:`~%s`'%attname for attname in self.extra_attributes]
                #rst +=', '.join()
            rst += self.format_list(attrs, indent=indent, columns=columns)
            rst += '\n'
        return rst

    def format_functions(self, indent=0, columns=None):
        """Format functions"""
        rst = ""
        if len(self.functions)+len(self.extra_functions):
            rst += self.indent(indent, ':Functions: ')
            funcs = []
            if self.functions:
                funcs += [self.format_ref(*func) for func in self.functions]
                #rst +=', '.join([self.format_ref(*func) for func in self.functions])
            if self.extra_functions:
                funcs += [':func:`~%s`'%funcname for funcname in self.extra_functions]
                #rst +=', '.join([':func:`~%s`'%funcname for funcname in self.extra_functions])
            rst += self.format_list(funcs, indent=indent, columns=columns)
            rst += '\n'
        return rst

    def format_classes(self, indent=0, columns=None):
        """Format classes

        .. todo:: Use :func:`inspect.getclasstree` or at least :func:`inspect.classify_class_attrs` in :class:`Overview`
        """

        rst = ""
        if len(self.classes)+len(self.extra_classes):
            rst += self.indent(indent, ':Classes: ')
            #rst += '\n'

            # Compact list or bullets?
#            if columns is None: columns = self.columns
#            nobj = max([len(])

            # Auto
            classes = []
            for clsname, cls in self.classes:
                kw = dict(clsname=clsname)
                #crst = self.indent(indent+1, self.format_ref(clsname, cls))
                crst = self.format_ref(clsname, cls)
                cls_attrs = self.class_contents[clsname]['attributes']
                cls_meths = self.class_contents[clsname]['methods']
                cls_xattrs = self.get_extra_class_attributes(clsname)
                cls_xmeths = self.get_extra_methods(clsname)
                if cls_attrs or cls_meths or cls_xattrs or cls_xmeths:
                    cls_attrs = [self.format_ref(*obj, **kw) for obj in cls_attrs]
                    cls_meths = [self.format_ref(*obj, **kw) for obj in cls_meths]
                    cls_xattrs = [':attr:`~%s`\n'%obj for obj in cls_xattrs]
                    cls_xmeths = [':meth:`~%s`\n'%obj for obj in cls_xmeths]
                    full = sorted(cls_attrs+cls_xattrs)+sorted(cls_meths+cls_xmeths)
                    crst += '\n'
                    crst += self.format_list(full, indent=indent+1, columns=columns, xindent='  ')
                    crst += '\n'
                classes.append(crst)

            # Extras
            for clsname in self.extra_classes:
                crst = self.indent(indent+1, ':class:`~%s`'%clsname)
                cls_xattrs = self.get_extra_class_attributes(clsname)
                cls_xmeths = self.get_extra_methods(clsname)
                if cls_xattrs or cls_xmeths:
                    cls_xattrs = [':attr:`~%s`\n'%obj for obj in cls_xattrs]
                    cls_xfuncs = [':func:`~%s`\n'%obj for obj in cls_xfuncs]
                    crst += '\n'
                    crst += self.format_list(sorted(cls_xattrs)+sorted(cls_xmeths),
                        indent=indent+1, columns=columns, xindent='  ')
                    crst += '\n'
                classes.append(crst)

            rst += self.format_list(classes, indent=indent, columns=1)

            rst += '\n'
        return rst

    def format(self, title_overview='Overview', title_content='Content', underline='-', indent=0, columns=None):
        """Format overview in rst format"""
        # Empty ?
        if not len(self.attributes+self.functions+self.classes+self.extra_attributes+
               self.extra_functions+self.extra_classes+self.extra_methods+self.extra_class_attributes):
            return ""
        rst = ""

        # Overview title
        rst += self.format_title(title_overview, underline, indent=indent)

        # Attributes
        rst += self.format_attributes(indent=indent, columns=columns)

        # Functions
        rst += self.format_functions(indent=indent, columns=columns)

        # Classes
        rst += self.format_classes(indent=indent, columns=columns)

        # Content title
        rst += self.format_title(title_content, underline, indent=indent)

        return rst




