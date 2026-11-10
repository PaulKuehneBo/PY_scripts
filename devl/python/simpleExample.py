import threading


def func(param, param2):
    """ func ist eine Test-Funktion 
    <param2>
        "bAddToPath":bool,
        "file":str
        "sID":int
    
    """
    d = dict()
    d.update( {"string": 88 })
    print("test {0}-{1}".format(param, param2))


class action:

    def __getattr__( self, name ):
        self.__dict__[name] = "neuer Eintrag"
    


func(4, 23)
func(42, 89)


verticies = (
    (1, -1, -1),
    (1, 1, -1),
    (-1, 1, -1),
    (-1, -1, -1),
    (1, -1, 1),
    (1, 1, 1),
    (-1, -1, 1),
    (-1, 1, 1)
)

edges = (
    (0, 1),
    (0, 3),
    (0, 4),
    (2, 1),
    (2, 3),
    (2, 7),
    (6, 3),
    (6, 4),
    (6, 7),
    (5, 1),
    (5, 4),
    (5, 7)
)




