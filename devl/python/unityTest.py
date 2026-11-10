##############################################################################################
##############################################################################################
from doctest import testfile
from unittest import result


G_TestRoutineDict = {}


# def registerTest(_funcTestMethod):
#     global G_TestRoutineDict
#     G_TestRoutineDict[_funcTestMethod] = []

#     def newDecorator(_TestSupplier):
#         def TestWrapper():
#             print(f"TestWrapper with {_TestSupplier}")
#             callableObject, dicArgs = _TestSupplier()

#             print(f"TestWrapper, now handling {callableObject} ( {dicArgs} )")
#             _funcTestMethod(callableObject, dicArgs)

#         G_TestRoutineDict[_funcTestMethod].append(TestWrapper)
#         return TestWrapper

#     return newDecorator


def registerTest(_funcTestMethod):
    global G_TestRoutineDict
    G_TestRoutineDict[_funcTestMethod] = []

    def newDecorator_withArg(_runsWithoutException, _showException=False):
        def newDecorator(_TestSupplier):
            def TestWrapper():
                print(f"TestWrapper with {_TestSupplier}")
                callableObject, dicArgs = _TestSupplier()

                exceptionThrown = False
                testFails = False
                try:
                    print(f"TestWrapper, now handling {callableObject} ( {dicArgs} )")
                    newDecorator_withArg.runsWithoutException = _runsWithoutException

                    _funcTestMethod(callableObject, dicArgs)
                except Exception as xEx:
                    exceptionThrown = True
                    if _runsWithoutException:
                        messageString = "Test Fails with Exception"
                        testFails = True
                    else:
                        messageString = (
                            "test runs sucessfully, exception raises as expected"
                        )

                if not exceptionThrown:
                    if _runsWithoutException:
                        messageString = "test runs sucessfully"
                    else:
                        messageString = "Test Fails, exception not thrown"
                        testFails = True

                print(messageString)
                return (not testFails, messageString)

            G_TestRoutineDict[_funcTestMethod].append(TestWrapper)
            # return TestWrapper

        return newDecorator

    return newDecorator_withArg


@registerTest
def testMethod(_callableObject, _dictArgs):

    print(f"testMethods calls {_callableObject} with {_dictArgs}")
    # calling with interface parameter,
    # checks the function and/or class c'tor,
    # it checks the runtime conditions
    tempValue = _callableObject(_dictArgs)

    # when not raising an exception, (otherwise the tempValue was not created)
    if testMethod.runsWithoutException:
        # the return value and/or the instance of the class
        print(f"creates tempObject {tempValue}, ready for checking")
        # when the function runs without exceptions,
        # a tempValue checking may causes exceptions
    else:
        # the exception was not raised, but expected
        # the tempValue can not be valid, a test doesn't make sense
        print(
            f"tempObject {tempValue}, is not valid,"
            " an exception was expected but not raised"
        )
        # rasing an exception yet, won't be necessary or correct
        # exceptions should only be thrown inside the callable


##############################################################################################
##############################################################################################
class CtestClass:
    def __init__(self, _dictArgs) -> None:
        self.d = _dictArgs


class CtestClass2:
    test: str

    def __init__(self, _dictArgs) -> None:
        self.dictionar = _dictArgs

    def m2( self, pa):
        pass

    @classmethod
    def c2( cls, pa):
        pass



@testMethod(False)
def class_ParameterProvider():
    print("42")
    CtestClass2.test = "test"
    dicArgs = {"test": 42}
    return (CtestClass, dicArgs)


@testMethod(True)
def class_ParameterProvider():
    print("3.1415")
    dicArgs = {"_Test_2_": 3.1415}
    return (CtestClass2, dicArgs)

class_ParameterProvider()

##############################################################################################
##############################################################################################
for d, functionList in G_TestRoutineDict.items():
    for r in functionList:
        result, messageString = r()
        if result is True:
            print(f"dieser Test reagiert wie erwartet: {messageString}")
        else:
            print(f"dieser Test reagiert nicht so wie erwartet: {messageString}")

##############################################################################################
##############################################################################################
