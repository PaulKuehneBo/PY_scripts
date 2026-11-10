from pyclbr import Class
import traceback

from collections import namedtuple
from typing import NamedTuple


class CCompletion(NamedTuple):
    a1: str
    a2: float


CField = NamedTuple("field", c1=str, c2=float)


A = CCompletion(8, 9)
F = CField("ee", "tf")
F.count


def AnnotatedTuple(typename, fields=None, /, **kwargs):
    class CTuple:
        _namedTuple = NamedTuple(typename, fields, **kwargs)

        def __init__(self, **kwargs) -> None:
            for key, val in kwargs.items():
                if isinstance(val, self._namedTuple.__annotations__[key]):
                    self.__setattr__(key, val)

    # dataclasses._set_new_attribute(cls, '__init__', dataclasses._create_fn(
    #     name='__init__', args=('self', 'dictArgs'),
    #     body=[
    #         'lFields_default = [xField for xField in fields(self) if xField.name.endswith("_default" )]',
    #         'lFields = [xField.name for xField in fields(self) if not xField.name.endswith("_default" )]',
    #         'for xField in lFields_default:',
    #         ' self.__setattr__(xField.name[:-8], xField.default)',
    #         'for key, val in dictArgs.items():', ' if key in lFields:',
    #         '  self.__setattr__(key, val)']))

    return CTuple


CAnnotatedTuple2 = NamedTuple("a2Tuple", buch=str, preis=int)
CAnnotatedTuple = AnnotatedTuple("annotationTuple", buch=str, preis=float)

karlMay = CAnnotatedTuple(buch="Winnetou", preis=19.99)
T2 = CAnnotatedTuple2(buch="wwww", preis=999)


class CAnnotatedTuple(NamedTuple):
    stueck: str
    preis: float


class Person:
    def __init__(self, _name, _age) -> None:
        self.name = _name
        self.age = _age


Dirk = Person("Dirk", 53)
Dirk_kind = Dirk
Dirk_kind.age = 12


line = namedtuple("line", ["stueck", "money"])
buch = line(money=25.95, stueck="Titel")
buch.stueck = "zwei"


def raiseEx():

    n = None
    n.func()
    return


try:
    raiseEx()
except Exception as Ex:
    print(f"das ging schief {Ex}")


class COwnEx(Exception):
    def __init__(self, sMsg) -> None:
        self.msg = sMsg
        pass


def raiseOwnEx():
    raise COwnEx("selber gemacht")


try:
    raiseOwnEx()
except Exception as Ex:
    print(f"Own Exception {Ex}")

print("Ende")
