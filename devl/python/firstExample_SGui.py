
import enum


class properties:
    _integer = 1
    _float   = 2
    _enum    = 3


class p:
    def __init__(self) -> None:
        self.type = "none"
        self.value = "none"
        self.hint = "none"


dim = p()
dim.type = "integer"
dim.value ="458"

prop = list()
prop.append( dim )

print( prop )


