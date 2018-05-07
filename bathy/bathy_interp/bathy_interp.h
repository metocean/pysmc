#include <Python.h>

PyObject * bathy_interpolate(PyObject * _file_name_in,
			     PyObject * _file_name_out,
			     PyObject * _lon0,
			     PyObject * _lon1,
			     PyObject * _dlon,
			     PyObject * _lat0,
			     PyObject * _lat1,
			     PyObject * _dlat,
			     PyObject * fill_value);
