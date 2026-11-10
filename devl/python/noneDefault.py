import enum
from  copy import deepcopy




def ifNoneDefault( f_arg, f_default):
    if f_arg is None:
        return deepcopy( f_default ) 
    
    return deepcopy(  f_arg )


class Bus:

    def __init__(self, passengers = None) -> None:
        self.passengers = ifNoneDefault( passengers, f_default=[])


    def entry( self, passengers = "Robert") -> None:
        if isinstance( passengers, str):
            self.passengers.append( passengers)
        elif isinstance( passengers, list):
            self.passengers += passengers

    def leave( self, passengers = "Chris") -> None:
        if isinstance( passengers, str):
            if passengers in self.passengers:
                self.passengers.remove( passengers)
        elif isinstance( passengers, list):
            self.passengers -= passengers            


alpha = enum.Enum( "Alpha", "alpha bety delta gamma")

class CDebug_simple:
    alpha = True

    def __init__(self) -> None:
        self.dirk = "dirk"
        pass


class CDebug(enum.Enum):
    alpha = True



def dPrint( f_msg ):
    if CDebug.alpha:
        print( f"debug Code {f_msg}\n")
    else:
        print( "release")


def function(inside=true):
    erstList =[2,42]
    d = Bus( erstList)
    i = 42

    dPrint( f"bus={Bus} = {i}")

    i = 89
    d.entry( passengers= ["2", "Robert"])
    dPrint( f"neue Ausgabe {Bus} = {i}")


function()


from functools import singledispatch


@singledispatch
def reactOnParam( param: object) -> str:
    return "default implementation"


@reactOnParam.register
def _(text: str)-> str:
    return f"string{text}"


@reactOnParam.register
def _(bus: Bus)-> str:
    return f"bus{bus.passengers}"


@reactOnParam.register
def _(f_dict: dict)-> str:
    return f"dictionary0{f_dict}"    


def func2():
    a = 42
    print( f"a: {reactOnParam(a)}\n" )
    a = "das ist ein String"
    print( f"str: {reactOnParam(a)}\n" )
    b = Bus( [42, "42-als-string"])
    print( f"Bus: {reactOnParam(b)}\n" )
    j = { "dirl":"rap", "claudia":"thies"}
    print( f"dict: {reactOnParam(j)}\n" )

    print( "Ende")    


func2()



from abc import get_cache_token
import types, weakref



##############################################################################################
##############################################################################################


WRAPPER_ASSIGNMENTS = ('__module__', '__name__', '__qualname__', '__doc__',
                       '__annotations__')
WRAPPER_UPDATES = ('__dict__',)
def update_wrapper(wrapper,
                   wrapped,
                   assigned = WRAPPER_ASSIGNMENTS,
                   updated = WRAPPER_UPDATES):
    """Update a wrapper function to look like the wrapped function

       wrapper is the function to be updated
       wrapped is the original function
       assigned is a tuple naming the attributes assigned directly
       from the wrapped function to the wrapper function (defaults to
       functools.WRAPPER_ASSIGNMENTS)
       updated is a tuple naming the attributes of the wrapper that
       are updated with the corresponding attribute from the wrapped
       function (defaults to functools.WRAPPER_UPDATES)
    """
    for attr in assigned:
        try:
            value = getattr(wrapped, attr)
        except AttributeError:
            pass
        else:
            setattr(wrapper, attr, value)
    for attr in updated:
        getattr(wrapper, attr).update(getattr(wrapped, attr, {}))
    # Issue #17482: set __wrapped__ last so we don't inadvertently copy it
    # from the wrapped function when updating __dict__
    wrapper.__wrapped__ = wrapped
    # Return the wrapper so this can be used as a decorator via partial()
    return wrapper


##############################################################################################
##############################################################################################


def _c3_merge(sequences):
    """Merges MROs in *sequences* to a single MRO using the C3 algorithm.

    Adapted from http://www.python.org/download/releases/2.3/mro/.

    """
    result = []
    while True:
        sequences = [s for s in sequences if s]   # purge empty sequences
        if not sequences:
            return result
        for s1 in sequences:   # find merge candidates among seq heads
            candidate = s1[0]
            for s2 in sequences:
                if candidate in s2[1:]:
                    candidate = None
                    break      # reject the current head, it appears later
            else:
                break
        if candidate is None:
            raise RuntimeError("Inconsistent hierarchy")
        result.append(candidate)
        # remove the chosen candidate
        for seq in sequences:
            if seq[0] == candidate:
                del seq[0]

def _c3_mro(cls, abcs=None):
    """Computes the method resolution order using extended C3 linearization.

    If no *abcs* are given, the algorithm works exactly like the built-in C3
    linearization used for method resolution.

    If given, *abcs* is a list of abstract base classes that should be inserted
    into the resulting MRO. Unrelated ABCs are ignored and don't end up in the
    result. The algorithm inserts ABCs where their functionality is introduced,
    i.e. issubclass(cls, abc) returns True for the class itself but returns
    False for all its direct base classes. Implicit ABCs for a given class
    (either registered or inferred from the presence of a special method like
    __len__) are inserted directly after the last ABC explicitly listed in the
    MRO of said class. If two implicit ABCs end up next to each other in the
    resulting MRO, their ordering depends on the order of types in *abcs*.

    """
    for i, base in enumerate(reversed(cls.__bases__)):
        if hasattr(base, '__abstractmethods__'):
            boundary = len(cls.__bases__) - i
            break   # Bases up to the last explicit ABC are considered first.
    else:
        boundary = 0
    abcs = list(abcs) if abcs else []
    explicit_bases = list(cls.__bases__[:boundary])
    abstract_bases = []
    other_bases = list(cls.__bases__[boundary:])
    for base in abcs:
        if issubclass(cls, base) and not any(
                issubclass(b, base) for b in cls.__bases__
            ):
            # If *cls* is the class that introduces behaviour described by
            # an ABC *base*, insert said ABC to its MRO.
            abstract_bases.append(base)
    for base in abstract_bases:
        abcs.remove(base)
    explicit_c3_mros = [_c3_mro(base, abcs=abcs) for base in explicit_bases]
    abstract_c3_mros = [_c3_mro(base, abcs=abcs) for base in abstract_bases]
    other_c3_mros = [_c3_mro(base, abcs=abcs) for base in other_bases]
    return _c3_merge(
        [[cls]] +
        explicit_c3_mros + abstract_c3_mros + other_c3_mros +
        [explicit_bases] + [abstract_bases] + [other_bases]
    )

def _compose_mro(cls, types):
    """Calculates the method resolution order for a given class *cls*.

    Includes relevant abstract base classes (with their respective bases) from
    the *types* iterable. Uses a modified C3 linearization algorithm.

    """
    bases = set(cls.__mro__)
    # Remove entries which are already present in the __mro__ or unrelated.
    def is_related(typ):
        return (typ not in bases and hasattr(typ, '__mro__')
                                 and issubclass(cls, typ))
    types = [n for n in types if is_related(n)]
    # Remove entries which are strict bases of other entries (they will end up
    # in the MRO anyway.
    def is_strict_base(typ):
        for other in types:
            if typ != other and typ in other.__mro__:
                return True
        return False
    types = [n for n in types if not is_strict_base(n)]
    # Subclasses of the ABCs in *types* which are also implemented by
    # *cls* can be used to stabilize ABC ordering.
    type_set = set(types)
    mro = []
    for typ in types:
        found = []
        for sub in typ.__subclasses__():
            if sub not in bases and issubclass(cls, sub):
                found.append([s for s in sub.__mro__ if s in type_set])
        if not found:
            mro.append(typ)
            continue
        # Favor subclasses with the biggest number of useful bases
        found.sort(key=len, reverse=True)
        for sub in found:
            for subcls in sub:
                if subcls not in mro:
                    mro.append(subcls)
    return _c3_mro(cls, abcs=mro)

def _find_impl(cls, registry):
    """Returns the best matching implementation from *registry* for type *cls*.

    Where there is no registered implementation for a specific type, its method
    resolution order is used to find a more generic implementation.

    Note: if *registry* does not contain an implementation for the base
    *object* type, this function may return None.

    """
    mro = _compose_mro(cls, registry.keys())
    match = None
    for t in mro:
        if match is not None:
            # If *match* is an implicit ABC but there is another unrelated,
            # equally matching implicit ABC, refuse the temptation to guess.
            if (t in registry and t not in cls.__mro__
                              and match not in cls.__mro__
                              and not issubclass(match, t)):
                raise RuntimeError("Ambiguous dispatch: {} or {}".format(
                    match, t))
            break
        if t in registry:
            match = t
    return registry.get(match)


##############################################################################################
##############################################################################################

def paramDispatch(func):
    """Single-dispatch generic function decorator.

    Transforms a function into a generic function, which can have different
    behaviours depending upon the type of its first argument. The decorated
    function acts as the default implementation, and additional
    implementations can be registered using the register() attribute of the
    generic function.
    """
    # There are many programs that use functools without singledispatch, so we
    # trade-off making singledispatch marginally slower for the benefit of
    # making start-up of such applications slightly faster.

    registry = {}
    dispatch_cache = dict()
    # dispatch_cache = weakref.WeakKeyDictionary()
    # cache_token = None

    def dispatch(cls):
        """generic_func.dispatch(cls) -> <function implementation>

        Runs the dispatch algorithm to return the best available implementation
        for the given *cls* registered on *generic_func*.

        """
        # nonlocal cache_token
        # if cache_token is not None:
        #     current_token = get_cache_token()
        #     if cache_token != current_token:
        #         dispatch_cache.clear()
        #         cache_token = current_token

        # try:
        #     impl = dispatch_cache[cls]
        # except KeyError:
        try:
            print( f"dispatch: {cls}")
            impl = registry[cls]
        except KeyError:
            impl = _find_impl(cls, registry)
        # dispatch_cache[cls] = impl
        return impl

    def register(cls, func=None):
        """generic_func.register(cls, func) -> func

        Registers a new implementation for the given *cls* on a *generic_func*.

        """
        # nonlocal cache_token
        if func is None:
            if isinstance(cls, type):
                return lambda f: register(cls, f)

            # only import typing if annotation parsing is necessary
            func = cls
            from typing import get_type_hints
            argClasses =  [ value_cls for varName, value_cls in get_type_hints(func).items() if varName is not "return"]
            if len( argClasses) != func.__code__.co_argcount:
                raise TypeError(
                    f"Invalid arguments to `register()`: {cls!r}. "
                    f"annotated args:{argClasses} differs from co_varnames:{func.__code__.co_varnames} "
                )

            argStr = "-".join( [ str(hash(c)) for c in argClasses])
            argTypes = [ not isinstance(cls,type) for cls in argClasses ]

            if any( argTypes ):
                falseTypes = [z[0] for z in zip( argClasses,argTypes) if z[1]]
                raise TypeError(
                    f"Invalid annotation for {falseTypes}. "
                )

        print( f"registration: {argStr}")
        registry[argStr] = func
        return func

    def wrapper(*args, **kw):
        if not args:
            raise TypeError(f'{funcname} requires at least '
                            '1 positional argument')
        argClasses = [ a.__class__ for a in args]
        argStr = "-".join( [ str(hash(c)) for c in argClasses])
        print( f"wrapper: {argStr}")
        # return dispatch(args[0].__class__)(*args, **kw)
        return dispatch(argStr)(*args, **kw)

    funcname = getattr(func, '__name__', 'singledispatch function')
    registry[object] = func
    wrapper.register = register
    wrapper.dispatch = dispatch
    wrapper.registry = types.MappingProxyType(registry)
    wrapper._clear_cache = dispatch_cache.clear
    update_wrapper(wrapper, func)
    return wrapper



##############################################################################################
##############################################################################################
def reactArg( *args ) -> str:
    return f"default implementation mit arg:{args}"
    pass

print( f"testObject: {reactArg( '980', 'abc', 'westen', '8765')}")

@paramDispatch
def react( *args ) -> str:
    return f"default implementation mit arg:{args}"


@react.register
def _( p1: str,p2: str ) -> str:
    return f"2 String-Parameter {p1}-{p2}"

@react.register
def _(text: str, i: int)-> str:
    return f"string_i{text} + {i}"


@react.register
def _(text_list: list, i: int, joinChar: str)-> str:
    return f"list { joinChar.join(text_list[:i])}"


@react.register
def _(f_dict: dict, f_key:str)-> str:
    return f"dictionary0{f_dict[f_key]}"    


def func3():
    print( f"testObject: {react( '980', 'abc', 'westen', '8765')}")
    print( f"int-int:   {react( '980', 8765)}")
    print( f"test-int:   {react( '980', 8765)}")
    j = { "dirl":"rap", "claudia":"thies"}
    print( f"test-dict:  {react( j, 'claudia')}")

func3()



print( "Ende")