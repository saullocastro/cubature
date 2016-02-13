#include "Python.h"

// Get the raw function pointer from a ctypes PyCFuncPtrObject (CFuncPtr)
void* get_ctypes_function_pointer(PyObject *obj);
