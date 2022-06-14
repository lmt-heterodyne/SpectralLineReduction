import numpy
from netCDF4 import Variable as NetCDFVariable
from lmtslr.utils.ordered_netcdf_dict import OrderedNetCDFDict, \
    OrderedHeaderDict
try:
    from collections import OrderedDict
except ImportError:
    OrderedDict = dict

class LMTHeader(OrderedNetCDFDict):
    """The most generic LMT Header class. This class is not usually
    instantiated by the user. It is part of the initialization of
    opening an LMT NetCDF file."""
    def __init__(self, ncvariables=None, dimensions=None,
                 suppress_items=['M1', 'Tiltmeter_0_', 'Tiltmeter_1_', 'XmlParser']):
        OrderedNetCDFDict.__init__(self) #, init_val=(), strict=False)
        #OrderedDict.__init__(self, init_val=(), strict=False)
        #self.ncvariables = ncvariables
        self.make_header_keys(ncvariables)
        self.dimensions = dimensions
        self.suppress_items = suppress_items

    def make_header_keys(self, ncvariables):
        heads = [name for name in ncvariables.keys() if name.find('Header') != -1]
        for head in heads:
            htype = head.split('.')[1]
            self[htype] = OrderedNetCDFDict()
            #self[htype] = OrderedDict()
        for key in self.keys():
            for subhead in [n for n in ncvariables.keys() if n.startswith('Header.%s.' % key)]:
                shead = subhead.split('.')[-1]
                self[key][shead] = ncvariables[subhead]

    def make_nominal_header(self):
        hdr = OrderedHeaderDict()
        for key, value in self.items():
            if key not in self.suppress_items:
                hdr[key] = OrderedHeaderDict()
                for k, v in value.items():
                    hdr[key][k] = self.get('%s.%s' % (key, k))
        return hdr                
    
    def _pprint_subhead(self, k, v):
        hstr = '{'
        for i, (key, val) in enumerate(v.items()):
            if isinstance(val, NetCDFVariable):
                if val.dtype == numpy.dtype('c'):
                    variab = val[:].tostring().strip()
                else:
                    variab = val[:]
                if i == 0:
                    hstr += "'%s': %s,\n" % (key, variab)
                else:
                    hstr += ' '*(len(k)+4)+"'%s': %s,\n" % (key, variab)
        hstr += '}\n'
        return hstr

    def utdate(self):
        """Returns python datetime equivalent for TimePlace.UTDate"""
        import datetime
        utd = self.get_scalar('TimePlace.UTDate')
        year = int(utd)
        if year % 4 == 0:
            days_in_year = 366.
        else:
            days_in_year = 365.
        secs = 24.*3600.*days_in_year*(utd-year)
        dt = datetime.datetime(year,1,1)+datetime.timedelta(seconds=secs)
        return dt



    
#     def __repr__(self):
#         hstr = ''
#         for k, v in self.items():
#             if type(v) == types.DictType:
#                 hstr += "{\n%s : %s\n}" % (k, self._pprint_subhead(k, v))
#         return hstr
