#include "Python.h"

// First part of PyCFuncPtrObject in Modules/_ctypes/ctypes.h
typedef struct {
    PyObject_HEAD char *b_ptr;
} cfuncptr_object;


// function pointer - any signature will work
typedef int (*_func_ptr)(double);

_func_ptr get_ctypes_function_pointer(PyObject *obj)
{
    return (*((void **) (((cfuncptr_object *)(obj))->b_ptr)));
}
