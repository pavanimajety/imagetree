# This file was automatically generated by SWIG (http://www.swig.org).
# Version 2.0.12
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.





from sys import version_info
if version_info >= (2,6,0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_Imagetree', [dirname(__file__)])
        except ImportError:
            import _Imagetree
            return _Imagetree
        if fp is not None:
            try:
                _mod = imp.load_module('_Imagetree', fp, pathname, description)
            finally:
                fp.close()
            return _mod
    _Imagetree = swig_import_helper()
    del swig_import_helper
else:
    import _Imagetree
del version_info
try:
    _swig_property = property
except NameError:
    pass # Python < 2.2 doesn't have 'property'.
def _swig_setattr_nondynamic(self,class_type,name,value,static=1):
    if (name == "thisown"): return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name,None)
    if method: return method(self,value)
    if (not static):
        self.__dict__[name] = value
    else:
        raise AttributeError("You cannot add attributes to %s" % self)

def _swig_setattr(self,class_type,name,value):
    return _swig_setattr_nondynamic(self,class_type,name,value,0)

def _swig_getattr(self,class_type,name):
    if (name == "thisown"): return self.this.own()
    method = class_type.__swig_getmethods__.get(name,None)
    if method: return method(self)
    raise AttributeError(name)

def _swig_repr(self):
    try: strthis = "proxy of " + self.this.__repr__()
    except: strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)

try:
    _object = object
    _newclass = 1
except AttributeError:
    class _object : pass
    _newclass = 0


NUMBER_OF_QUADRANTS = _Imagetree.NUMBER_OF_QUADRANTS
NorthEast = _Imagetree.NorthEast
SouthEast = _Imagetree.SouthEast
SouthWest = _Imagetree.SouthWest
NorthWest = _Imagetree.NorthWest
class Pixel(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Pixel, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Pixel, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _Imagetree.new_Pixel(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_setmethods__["R"] = _Imagetree.Pixel_R_set
    __swig_getmethods__["R"] = _Imagetree.Pixel_R_get
    if _newclass:R = _swig_property(_Imagetree.Pixel_R_get, _Imagetree.Pixel_R_set)
    __swig_setmethods__["G"] = _Imagetree.Pixel_G_set
    __swig_getmethods__["G"] = _Imagetree.Pixel_G_get
    if _newclass:G = _swig_property(_Imagetree.Pixel_G_get, _Imagetree.Pixel_G_set)
    __swig_setmethods__["B"] = _Imagetree.Pixel_B_set
    __swig_getmethods__["B"] = _Imagetree.Pixel_B_get
    if _newclass:B = _swig_property(_Imagetree.Pixel_B_get, _Imagetree.Pixel_B_set)
    __swig_setmethods__["intensity"] = _Imagetree.Pixel_intensity_set
    __swig_getmethods__["intensity"] = _Imagetree.Pixel_intensity_get
    if _newclass:intensity = _swig_property(_Imagetree.Pixel_intensity_get, _Imagetree.Pixel_intensity_set)
    __swig_destroy__ = _Imagetree.delete_Pixel
    __del__ = lambda self : None;
Pixel_swigregister = _Imagetree.Pixel_swigregister
Pixel_swigregister(Pixel)

class Imagetree(_object):
    __swig_setmethods__ = {}
    __setattr__ = lambda self, name, value: _swig_setattr(self, Imagetree, name, value)
    __swig_getmethods__ = {}
    __getattr__ = lambda self, name: _swig_getattr(self, Imagetree, name)
    def __init__(self, *args, **kwargs): raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _Imagetree.delete_Imagetree
    __del__ = lambda self : None;
    def isLeaf(self): return _Imagetree.Imagetree_isLeaf(self)
    def isNode(self): return _Imagetree.Imagetree_isNode(self)
    def numberOfLeaves(self): return _Imagetree.Imagetree_numberOfLeaves(self)
    def numberOfNodes(self): return _Imagetree.Imagetree_numberOfNodes(self)
    def numberOfSubTrees(self): return _Imagetree.Imagetree_numberOfSubTrees(self)
    def value(self): return _Imagetree.Imagetree_value(self)
Imagetree_swigregister = _Imagetree.Imagetree_swigregister
Imagetree_swigregister(Imagetree)

class QuadLeaf(Imagetree):
    __swig_setmethods__ = {}
    for _s in [Imagetree]: __swig_setmethods__.update(getattr(_s,'__swig_setmethods__',{}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, QuadLeaf, name, value)
    __swig_getmethods__ = {}
    for _s in [Imagetree]: __swig_getmethods__.update(getattr(_s,'__swig_getmethods__',{}))
    __getattr__ = lambda self, name: _swig_getattr(self, QuadLeaf, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _Imagetree.new_QuadLeaf(*args)
        try: self.this.append(this)
        except: self.this = this
    def isLeaf(self): return _Imagetree.QuadLeaf_isLeaf(self)
    def numberOfNodes(self): return _Imagetree.QuadLeaf_numberOfNodes(self)
    def numberOfLeaves(self): return _Imagetree.QuadLeaf_numberOfLeaves(self)
    def numberOfSubTrees(self): return _Imagetree.QuadLeaf_numberOfSubTrees(self)
    def value(self): return _Imagetree.QuadLeaf_value(self)
    def son(self, *args): return _Imagetree.QuadLeaf_son(self, *args)
    __swig_destroy__ = _Imagetree.delete_QuadLeaf
    __del__ = lambda self : None;
QuadLeaf_swigregister = _Imagetree.QuadLeaf_swigregister
QuadLeaf_swigregister(QuadLeaf)

class QuadNode(Imagetree):
    __swig_setmethods__ = {}
    for _s in [Imagetree]: __swig_setmethods__.update(getattr(_s,'__swig_setmethods__',{}))
    __setattr__ = lambda self, name, value: _swig_setattr(self, QuadNode, name, value)
    __swig_getmethods__ = {}
    for _s in [Imagetree]: __swig_getmethods__.update(getattr(_s,'__swig_getmethods__',{}))
    __getattr__ = lambda self, name: _swig_getattr(self, QuadNode, name)
    __repr__ = _swig_repr
    def __init__(self, *args): 
        this = _Imagetree.new_QuadNode(*args)
        try: self.this.append(this)
        except: self.this = this
    __swig_destroy__ = _Imagetree.delete_QuadNode
    __del__ = lambda self : None;
    def isLeaf(self): return _Imagetree.QuadNode_isLeaf(self)
    def numberOfLeaves(self): return _Imagetree.QuadNode_numberOfLeaves(self)
    def numberOfNodes(self): return _Imagetree.QuadNode_numberOfNodes(self)
    def numberOfSubTrees(self): return _Imagetree.QuadNode_numberOfSubTrees(self)
    def value(self): return _Imagetree.QuadNode_value(self)
    def son(self, *args): return _Imagetree.QuadNode_son(self, *args)
    def killSon(self, *args): return _Imagetree.QuadNode_killSon(self, *args)
QuadNode_swigregister = _Imagetree.QuadNode_swigregister
QuadNode_swigregister(QuadNode)

# This file is compatible with both classic and new-style classes.


