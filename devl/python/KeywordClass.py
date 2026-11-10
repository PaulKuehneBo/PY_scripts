#!/usr/bin/env python3
# -*- coding:utf-8 -*-
###
# File: \CheckFuncTypes.py
# Created Date: Wednesday, May 25th 2022, 11:51:33 am
# Author: Christian Perwass (CR/AEC5)
# <LICENSE id="Apache-2.0">
#
#   Image-Render Base Functions module
#   Copyright 2022 Robert Bosch GmbH and its subsidiaries
#
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.
#
# </LICENSE>
###

from multiprocessing.spawn import old_main_modules
import sys
from dataclasses import dataclass, fields
import dataclasses
from xml.etree.ElementInclude import default_loader


def _processKeywordClass(cls):
    """this wrapper is similiar to dataclasses implementation of the python lib"""
    # Now that dicts retain insertion order, there's no reason to use
    # an ordered dict.  I am leveraging that ordering here, because
    # derived class fields overwrite base class fields, but the order
    # is defined by the base class, which is found first.
    fields = {}

    if cls.__module__ in sys.modules:
        globals = sys.modules[cls.__module__].__dict__
    else:
        # Theoretically this can happen if someone writes
        # a custom string to cls.__module__.  In which case
        # such dataclass won't be fully introspectable
        # (w.r.t. typing.get_type_hints) but will still function
        # correctly.
        globals = {}

    cls_annotations = cls.__dict__.get('__annotations__', {})

    # Now find fields in our class.  While doing so, validate some
    # things, and set the default values (as class attributes) where
    # we can.
    cls_fields = [dataclasses._get_field(cls, name, type)
                  for name, type in cls_annotations.items()]
    for f in cls_fields:
        name = f.name

        fields[name] = f
        fields[name].default = name

        setattr(cls, name, name)

    # Do we have any Field members that don't also have annotations?
    for name, value in cls.__dict__.items():
        if isinstance(value, dataclasses.Field) and not name in cls_annotations:
            raise TypeError(f'{name!r} is a field but has no type annotation')

    # Remember all of the fields on our class (including bases).  This
    # also marks this class as being a dataclass.
    setattr(cls, dataclasses._FIELDS, fields)

    dataclasses._set_new_attribute(cls, '__init__', dataclasses._create_fn(
        name='__init__', args=('self', 'dictArgs'),
        body=[
            'lFields_default = [xField for xField in fields(self) if xField.name.endswith("_default" )]',
            'lFields = [xField.name for xField in fields(self) if not xField.name.endswith("_default" )]',
            'for xField in lFields_default:',
            ' self.__setattr__(xField.name[:-8], xField.default)',
            'for key, val in dictArgs.items():', ' if key in lFields:',
            '  self.__setattr__(key, val)']))

    return cls


def keywordclass(cls=None):
    """
    this wrapper is similiar to dataclasses implementation of the python lib

    Returns the same class as was passed in, with dunder methods
    added based on the fields defined in the class.
    """

    def wrap(cls):
        return _processKeywordClass(cls)

    # See if we're being called as @keywordclass or @keywordclass().
    if cls is None:
        # We're called with parens.
        return wrap

    # We're called as @keywordclass without parens.
    return wrap(cls)


class CReadOnly(type):
    def __setattr__(self, __name: str, __value) -> None:
        raise ValueError(f"{self}: do not change {__name}")


def _processConstClass(cls):
    """this wrapper is similiar to dataclasses implementation of the python lib"""
    # Now that dicts retain insertion order, there's no reason to use
    # an ordered dict.  I am leveraging that ordering here, because
    # derived class fields overwrite base class fields, but the order
    # is defined by the base class, which is found first.

    fields = {}

    if cls.__module__ in sys.modules:
        globals = sys.modules[cls.__module__].__dict__
    else:
        # Theoretically this can happen if someone writes
        # a custom string to cls.__module__.  In which case
        # such dataclass won't be fully introspectable
        # (w.r.t. typing.get_type_hints) but will still function
        # correctly.
        globals = {}

    cls_annotations = cls.__dict__.get('__annotations__', {})

    # Now find fields in our class.  While doing so, validate some
    # things, and set the default values (as class attributes) where
    # we can.
    cls_fields = [dataclasses._get_field(cls, name, type)
                  for name, type in cls_annotations.items()]
    for f in cls_fields:
        name = f.name

        fields[name] = f
        fields[name].default = name

        setattr(cls, name, name)

    # Do we have any Field members that don't also have annotations?
    for name, value in cls.__dict__.items():
        if isinstance(value, dataclasses.Field) and not name in cls_annotations:
            raise TypeError(f'{name!r} is a field but has no type annotation')

    # Remember all of the fields on our class (including bases).  This
    # also marks this class as being a dataclass.
    setattr(cls, dataclasses._FIELDS, fields)

    # dataclasses._set_new_attribute(cls, '__init__', dataclasses._create_fn(
    #     name='__init__', args=('self', 'dictArgs'),
    #     body=[
    #         'lFields_default = [xField for xField in fields(self) if xField.name.endswith("_default" )]',
    #         'lFields = [xField.name for xField in fields(self) if not xField.name.endswith("_default" )]',
    #         'for xField in lFields_default:',
    #         ' self.__setattr__(xField.name[:-8], xField.default)',
    #         'for key, val in dictArgs.items():', ' if key in lFields:',
    #         '  self.__setattr__(key, val)']))
    class CNewClass(cls, metaclass=CReadOnly):
        pass

    return CNewClass


def constclass(cls=None):
    """
    this wrapper is similiar to dataclasses implementation of the python lib

    Returns the same class as was passed in, with dunder methods
    added based on the fields defined in the class.
    """

    def wrap(cls):
        return _processConstClass(cls)

    # See if we're being called as @keywordclass or @keywordclass().
    if cls is None:
        # We're called with parens.
        return wrap

    # We're called as @keywordclass without parens.
    return wrap(cls)


@keywordclass
class CKeyWordParam:  # (CParamAssign):
    sDTI: str
    xModule: str = None
    sUnit: str = "rad"
    sVeryImportant = None

    def GetModule(self):
        return

    pass


@constclass
class CKeyWordParam2:

    sDTI: str
    xModule: str = None
    sUnit: str = "rad"
    sVeryImportant = None

    def GetModule(cls):
        return 6

    pass


@constclass
class CKWParam1:

    sDTI: str
    xBlender: str = None
    fPixel: str = "rad"
    sVeryImportant = None

    pass


@constclass
class CKWParam2:

    sXML: str
    xBase: str = None
    sLense: str = "rad"
    sVeryImportant = None

    pass


def getValue(f_dictArgs: dict, f_xKey, f_default=None):
    if isinstance(f_xKey, str):
        return f_dictArgs.get(f_xKey, f_default)

    sKey = str(f_xKey)
    return f_dictArgs.get(sKey, f_default)


def makeInterface(xFunc, lField):
    for sField in lField:
        setattr(xFunc, sField, sField)


def DeltaRotationEuler(dicArgs):
    """test"""
    print("test")

    data = CParam(dicArgs)

    print(f"{data.sUnit}")
    print(f"{ DeltaRotationEuler.sDTI} ")
    pass


makeInterface(DeltaRotationEuler, ["sDTI", "xValue", "sUnit"])


DeltaRotationEuler.interface = {
    "value": "blender/animate/object/rotate/const:1",
    "type": "label",
    "hint": "entry point identification",
    "sUnit": "entry point identification",
}


if __name__ == "__main__":
    dicArgs = {
        "value": "blender/animate/object/rotate/const:1",
        "type": "label",
        "hint": "entry point identification",
        "sUnit": "entry point identification",
        "sVeryImportant": "wird das wirklich überschrieben",
    }

    DeltaRotationEuler(dicArgs=dicArgs)
