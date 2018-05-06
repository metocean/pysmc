#include <Python.h>


//Actual module method definition - this is the code that will be called by
//hello_module.print_hello_world
PyObject * bathy_interp_interface       (PyObject *self,
					 PyObject *args,
					 PyObject *keywds);


//Method definition object for this extension, these argumens mean:
//ml_name: The name of the method
//ml_meth: Function pointer to the method implementation
//ml_flags: Flags indicating special features of this method, such as
//          accepting arguments, accepting keyword arguments, being a
//          class method, or being a static method of a class.
//ml_doc:  Contents of this method's docstring
static PyMethodDef bathy_interp_module_methods[] = { 
  {   
    "bathy_interp",
    bathy_interp_interface,
    METH_VARARGS| METH_KEYWORDS,
    "Interpolates a bathymetry."
  },
  {NULL, NULL, 0, NULL}
};

//Module definition
//The arguments of this structure tell Python what to call your extension,
//what it's methods are and where to look for it's method definitions
static struct PyModuleDef bathy_interp_module_definition = { 
    PyModuleDef_HEAD_INIT,
    "bathy_interp",
    "A Python module that does bathymetry interpolation",
    -1, 
    bathy_interp_module_methods
};

//Module initialization
//Python calls this function when importing your extension. It is important
//that this function is named PyInit_[[your_module_name]] exactly, and matches
//the name keyword argument in setup.py's setup() call.
PyMODINIT_FUNC PyInit_bathy_interp(void)
{
    Py_Initialize();

    return PyModule_Create(&bathy_interp_module_definition);
}
