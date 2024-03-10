#define PY_SSIZE_T_CLEAN
#define _USE_MATH_DEFINES
#include <stdio.h>
#include <stdlib.h>
#include <Python.h>
#include <Windows.h>
#include <math.h>
#ifndef M_PI
#define M_PI (3.14159265358979323846264338327950288)
#endif


long long fib(int n)
{
    if (n == 0 || n == 1)
        return 1;
    else if (n > 1)
        return fib(n - 1) + fib(n - 2);
    else
    {
        fprintf(stderr, "Error! Argument should be non-negative number.\n");
        exit(EXIT_FAILURE);
    }
}


double normal_density(double x)
{
    return 1 / sqrt(2 * M_PI) * exp(-1 * pow(x, 2) / 2);
}


double normal_density_gen(double x, double std, double mean)
{
    if (std <= 0)
    {
        fprintf(stderr, "Error! Std must be greater than zero.\n");
        exit(EXIT_FAILURE);
    }
    return 1 / sqrt(2 * M_PI * std) * exp(-1 * pow(x - mean, 2) / (2 * pow(std, 2)));
}


double uniform_density(double x, double a, double b)
{
    if (a >= b)
    {
        fprintf(stderr, "Error! Parameter a must be smaller than b.\n");
        exit(EXIT_FAILURE);
    }
    if (a <= x && x <= b)
        return 1 / (b - a);
    else
        return 0;
}


double tstudent_density(double x, int degree_of_freedom)
{
    if (degree_of_freedom < 0)
    {
        fprintf(stderr, "Error! Argument 'degree_of_freedom' should be non-negative number.\n");
        exit(EXIT_FAILURE);
    }

    return tgamma((degree_of_freedom + 1) / 2) / (sqrt(M_PI * degree_of_freedom) * tgamma(degree_of_freedom / 2)) * pow(1 + pow(x, 2) / degree_of_freedom, -0.5 * (degree_of_freedom + 1));
}


double exp_density(double x, double lambda)
{
    if (lambda < 0)
    {
        fprintf(stderr, "Error! Argument 'lambda' should be positive number.\n");
        exit(EXIT_FAILURE);
    }

    if (x <= 0)
        return 0;
    else
        return lambda * exp(-lambda * x);
}


double Cauchy_density(double x, double shift, double gamma)
{
    if (gamma < 0)
    {
        fprintf(stderr, "Error! Argument 'gamma' should be positive number.\n");
        exit(EXIT_FAILURE);
    }
    return gamma / (M_PI * (pow(x - shift, 2) + pow(gamma, 2)));
}


double chi_squared_density(double x, int degree_of_freedom)
{
    if (degree_of_freedom < 0)
    {
        fprintf(stderr, "Error! Argument 'degree_of_freedom' should be non-negative number.\n");
        exit(EXIT_FAILURE);
    }

    if (x > 0)
        return pow(x, degree_of_freedom - 1) * exp(-0.5 * pow(x, 2)) / (pow(2, degree_of_freedom / 2 - 1) * tgamma(degree_of_freedom / 2));
    else
        return 0.0;
}


struct cdf_data
{
    double* args;
    double* probs;
    double cdf;
    int size;
};


struct cdf_data  normal_cdf(double x)
{
    struct cdf_data  fin;
    int num_steps;
    if (x > -3.0)
        num_steps = (int)((x + 3.0) * 100);
    else
        num_steps = 1;
    fin.args = malloc(sizeof(double) * num_steps);
    fin.probs = malloc(sizeof(double) * num_steps);
    if (x > -3.0)
    {
        int iterator = 0;
        for (double i = -3.0; i <= x; i += 0.01)
        {
            fin.args[iterator] = i;
            if (iterator > 0)
                fin.probs[iterator] = fin.probs[iterator - 1] + normal_density(i) * 0.01;
            else
                fin.probs[iterator] = normal_density(i) * 0.01;
            iterator++;
        }
        fin.cdf = fin.probs[iterator - 1];
        fin.size = iterator;
        return fin;
    }
    else
    {
        fin.args[0] = x;
        fin.probs[0] = 0.0;
        fin.cdf = 0.0;
        fin.size = 1;
        return fin;
    }
}


struct cdf_data  normal_gen_cdf(double x, double mean, double std)
{
    double min_limit = -3.0 + mean;
    struct cdf_data  fin;
    int num_steps;
    if (x > min_limit)
        num_steps = (int)((x - min_limit) * 100);
    else
        num_steps = 1;
    fin.args = malloc(sizeof(double) * num_steps);
    fin.probs = malloc(sizeof(double) * num_steps);
    if (x > min_limit)
    {
        int iterator = 0;
        for (double i = min_limit; i <= x; i += 0.01)
        {
            fin.args[iterator] = i;
            if (iterator > 0)
                fin.probs[iterator] = fin.probs[iterator - 1] + normal_density_gen(i, std, mean) * 0.01;
            else
                fin.probs[iterator] = normal_density_gen(i, std, mean) * 0.01;
            iterator++;
        }
        fin.cdf = fin.probs[iterator - 1];
        fin.size = iterator;
        return fin;
    }
    else
    {
        fin.args[0] = x;
        fin.probs[0] = 0.0;
        fin.cdf = 0.0;
        fin.size = 1;
        return fin;
    }
}


struct cdf_data  uniform_cdf(double x, double a, double b)
{
    if (a >= b)
    {
        fprintf(stderr, "Error! Parameter a must be smaller than b.\n");
        exit(EXIT_FAILURE);
    }
    struct cdf_data  fin;
    int num_steps = (int)((x + a) * 100);
    fin.args = malloc(sizeof(double) * num_steps);
    fin.probs = malloc(sizeof(double) * num_steps);

    if (x > a)
    {
        int iterator = 0;
        for (double i = a; i <= x; i += 0.01)
        {
            fin.args[iterator] = i;
            if (iterator > 0)
                fin.probs[iterator] = (i + 0.01 - a) / (b - a);
            else
                fin.probs[iterator] = 0;
            iterator++;
        }
        fin.cdf = fin.probs[iterator - 1];
        fin.size = iterator;
        return fin;
    }
    else
    {
        fin.args[0] = x;
        fin.probs[0] = 0.0;
        fin.cdf = 0.0;
        fin.size = 1;
        return fin;
    }
}


struct cdf_data  tstudent_cdf(double x, int degree_of_freedom)
{
    if (degree_of_freedom < 0)
    {
        fprintf(stderr, "Error! Argument 'degree_of_freedom' should be non-negative number.\n");
        exit(EXIT_FAILURE);
    }
    struct cdf_data  fin;
    int num_steps = (int)((x + 10) * 100);
    fin.args = malloc(sizeof(double) * num_steps);
    fin.probs = malloc(sizeof(double) * num_steps);
    int iterator = 0;
    for (double i = -10.0; i <= x; i += 0.01)
    {
        fin.args[iterator] = i;
        if (iterator > 0)
            fin.probs[iterator] = fin.probs[iterator - 1] + tstudent_density(i, degree_of_freedom) * 0.01;
        else
            fin.probs[iterator] = tstudent_density(i, degree_of_freedom) * 0.01;
        iterator++;
    }
    fin.cdf = fin.probs[iterator - 1];
    fin.size = iterator;
    return fin;
}


struct cdf_data  exp_cdf(double x, double lambda)
{
    if (lambda < 0)
    {
        fprintf(stderr, "Error! Argument 'lambda' should be positive number.\n");
        exit(EXIT_FAILURE);
    }
    struct cdf_data  fin;
    int num_steps;
    if (x > 0)
        num_steps = (int)(x * 100);
    else
        num_steps = 1;
    fin.args = malloc(sizeof(double) * num_steps);
    fin.probs = malloc(sizeof(double) * num_steps);
    if (x > 0)
    {
        int iterator = 0;
        for (double i = 0.0; i <= x; i += 0.01)
        {
            fin.args[iterator] = i;
            if (iterator > 0)
                fin.probs[iterator] = fin.probs[iterator - 1] + (1 - exp(-1 * lambda * i)) * 0.01;
            else
                fin.probs[iterator] = (1 - exp(-1 * lambda * i)) * 0.01;
            iterator++;
        }
        fin.cdf = fin.probs[iterator - 1];
        fin.size = iterator;
        return fin;
    }
    else
    {
        fin.args[0] = x;
        fin.probs[0] = 0.0;
        fin.cdf = 0.0;
        fin.size = 1;
        return fin;
    }
}


struct cdf_data  Cauchy_cdf(double x, double shift, double gamma)
{
    if (gamma < 0)
    {
        fprintf(stderr, "Error! Argument 'gamma' should be positive number.\n");
        exit(EXIT_FAILURE);
    }
    double min_value = -20.0 + shift;
    struct cdf_data  fin;
    int num_steps;
    if (x > min_value)
        num_steps = (int)((x - min_value) * 100);
    else
        num_steps = 1;
    fin.args = malloc(sizeof(double) * num_steps);
    fin.probs = malloc(sizeof(double) * num_steps);
    int iterator = 0;
    if (x > min_value)
    {
        for (double i = min_value; i <= x; i += 0.01)
        {
            fin.args[iterator] = i;
            if (iterator > 0)
                fin.probs[iterator] = fin.probs[iterator - 1] + (0.5 + 1 / M_PI * atan((i - shift) / gamma)) * 0.01;
            else
                fin.probs[iterator] = (0.5 + 1 / M_PI * atan((i - shift) / gamma)) * 0.01;
            iterator++;
        }
        fin.cdf = fin.probs[iterator - 1];
        fin.size = iterator;
        return fin;
    }
    else
    {
        fin.args[0] = x;
        fin.probs[0] = 0.0;
        fin.cdf = 0.0;
        fin.size = 1;
        return fin;
    }
}


struct cdf_data  chi_squared_cdf(double x, int degree_of_freedom)
{
    if (degree_of_freedom < 0)
    {
        fprintf(stderr, "Error! Argument 'degree_of_freedom' should be non-negative number.\n");
        exit(EXIT_FAILURE);
    }
    struct cdf_data  fin;
    int num_steps;
    if (x > 0)
        num_steps = (int)(x * 100);
    else
        num_steps = 1;
    fin.args = malloc(sizeof(double) * num_steps);
    fin.probs = malloc(sizeof(double) * num_steps);
    if (x > 0)
    {
        int iterator = 0;
        for (double i = 0.0; i <= x; i += 0.01)
        {
            fin.args[iterator] = i;
            if (iterator > 0)
                fin.probs[iterator] = fin.probs[iterator - 1] + chi_squared_density(i, degree_of_freedom) * 0.01;
            else
                fin.probs[iterator] = chi_squared_density(i, degree_of_freedom) * 0.01;
            iterator++;
        }
        fin.cdf = fin.probs[iterator - 1];
        fin.size = iterator;
        return fin;
    }
    else
    {
        fin.args[0] = x;
        fin.probs[0] = 0.0;
        fin.cdf = 0.0;
        fin.size = 1;
        return fin;
    }
}



double find_percentile(double x, double prob, double arg_minimal, struct cdf_data  norm)
{
    double arg_min = arg_minimal;
    double dist = 1.0;
    int iterator = 0;

    for (double i = arg_minimal; i < x; i += 0.01)
    {
        if (fabs(norm.probs[iterator] - prob) < dist)
        {
            dist = fabs(norm.probs[iterator] - prob);
            arg_min = norm.args[iterator];
        }
        iterator++;
    }
    free(norm.args);
    free(norm.probs);

    return arg_min;
}



double percentile_normal_cdf(double prob)
{
    if (prob == 0)
        return -3.0;

    double x = -2.0;
    struct cdf_data  norm = normal_cdf(x);

    while (norm.cdf < prob)
    {
        x += 0.5;
        norm = normal_cdf(x);
    }
    return find_percentile(x, prob, -3.0, norm);
}


double percentile_normal_gen_cdf(double prob, double mean, double std)
{
    double min_limit = -3.0 + mean;
    if (prob == 0)
        return -3.0 + mean;

    double x = min_limit + 1.0;
    struct cdf_data  norm = normal_gen_cdf(x, mean, std);

    while (norm.cdf < prob)
    {
        x += 0.5;
        norm = normal_gen_cdf(x, mean, std);
    }

    return find_percentile(x, prob, min_limit, norm);
}


double percentile_uniform_cdf(double prob, double a, double b)
{

    if (prob == 0)
        return a;

    double x = a + (b - a) / 4;
    struct cdf_data  norm = uniform_cdf(x, a, b);

    while (norm.cdf < prob)
    {
        x += (b - a) / 4;
        norm = uniform_cdf(x, a, b);
    }

    return find_percentile(x, prob, a, norm);
}


double percentile_tstudent_cdf(double prob, int degree_of_freedom)
{

    if (prob == 0)
        return -10.0;

    double x = -9.0;
    struct cdf_data norm = tstudent_cdf(x, degree_of_freedom);

    while (norm.cdf < prob)
    {
        x += 0.5;
        norm = tstudent_cdf(x, degree_of_freedom);
    }

    return find_percentile(x, prob, -10.0, norm);
}


double percentile_exp_cdf(double prob, double lambda)
{

    if (prob == 0)
        return 0.0;

    double x = 1.0;
    struct cdf_data  norm = exp_cdf(x, lambda);

    while (norm.cdf < prob)
    {
        x += 0.5;
        norm = exp_cdf(x, lambda);
    }

    return find_percentile(x, prob, 0.0, norm);
}


double percentile_Cauchy_cdf(double prob, double shift, double gamma)
{

    double min_value = -20.0 + shift;
    if (prob == 0)
        return min_value;

    double x = -19.0 + shift;
    struct cdf_data  norm = Cauchy_cdf(x, shift, gamma);

    while (norm.cdf < prob)
    {
        x += 0.5;
        norm = Cauchy_cdf(x, shift, gamma);
    }

    return find_percentile(x, prob, min_value, norm);
}


double percentile_chi_squared_cdf(double prob, int degree_of_freedom)
{

    if (prob == 0)
        return 0.0;

    double x = 1.0;
    struct cdf_data  norm = chi_squared_cdf(x, degree_of_freedom);

    while (norm.cdf < prob)
    {
        x += 0.5;
        norm = chi_squared_cdf(x, degree_of_freedom);
    }

    return find_percentile(x, prob, 0.0, norm);
}



static PyObject* normal_density_py(PyObject* self, PyObject* args)
{
    PyObject* result = NULL;
    double x;
    if (!PyArg_ParseTuple(args, "d", &x))
        return NULL;
    return Py_BuildValue("d", normal_density(x));
}



static PyObject* density_py(PyObject* self, PyObject* args, PyObject* kwargs)
{
    PyObject* result = NULL;
    double x;
    char* density_type = NULL;
    double parameter_1 = 0.0;
    double parameter_2 = 0.0;
    static char* keywords[] = { "x", "density_type", "parameter_1", "parameter_2", NULL };

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ds|dd", &keywords, &x, &density_type, &parameter_1, &parameter_2))
        return NULL;
    //printf("\n x: %f\n density type: %s\n parameter_1: %f\n parameter_2: %f\n parameter_3: %d\n", x, density_type, parameter_1, parameter_2, parameter_3);
    if (strcmp(density_type, "normal") != 0 && strcmp(density_type, "cauchy") != 0 && strcmp(density_type, "uniform") != 0 &&
        strcmp(density_type, "tstudent") != 0 && strcmp(density_type, "chi_squared") != 0 && strcmp(density_type, "exp") != 0)
    {
        PyErr_SetString(PyExc_ValueError, "Unknown density type");
        return NULL;
    }

    if (strcmp(density_type, "normal") == 0)
        return Py_BuildValue("d", normal_density_gen(x, parameter_2, parameter_1));
    else if (strcmp(density_type, "uniform") == 0)
        return Py_BuildValue("d", uniform_density(x, parameter_1, parameter_2));
    else if (strcmp(density_type, "tstudent") == 0)
        return Py_BuildValue("d", tstudent_density(x, (int)parameter_1));
    else if (strcmp(density_type, "chi_squared") == 0)
        return Py_BuildValue("d", chi_squared_density(x, (int)parameter_1));
    else if (strcmp(density_type, "exp") == 0)
        return Py_BuildValue("d", exp_density(x, parameter_1));
    else if (strcmp(density_type, "cauchy") == 0)
        return Py_BuildValue("d", Cauchy_density(x, parameter_1, parameter_2));
    else
    {
        fprintf(stderr, "Error! Allowed types of density are: ['normal', 'uniform', 'tstudent', 'chi_squared', 'exp', 'cauchy']!\n");
        exit(EXIT_FAILURE);
    }
}


static PyObject* cdf_py(PyObject* self, PyObject* args, PyObject* kwargs)
{
    PyObject* result = NULL;
    double x;
    char* cdf_type = NULL;
    double parameter_1 = 0.0;
    double parameter_2 = 0.0;
    struct cdf_data  norm;
    PyObject* dict_res = PyDict_New();
    PyObject* args_list = PyList_New(0);
    PyObject* probs_list = PyList_New(0);

    static char* keywords[] = { "x", "cdf_type", "parameter_1", "parameter_2", NULL };

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ds|dd", &keywords, &x, &cdf_type, &parameter_1, &parameter_2))
        return NULL;
    //printf("\n x: %f\n density type: %s\n parameter_1: %f\n parameter_2: %f\n parameter_3: %d\n", x, density_type, parameter_1, parameter_2, parameter_3);
    if (strcmp(cdf_type, "normal") != 0 && strcmp(cdf_type, "cauchy") != 0 && strcmp(cdf_type, "uniform") != 0 &&
        strcmp(cdf_type, "tstudent") != 0 && strcmp(cdf_type, "chi_squared") != 0 && strcmp(cdf_type, "exp") != 0)
    {
        PyErr_SetString(PyExc_ValueError, "Unknown density type");
        return NULL;
    }

    if (strcmp(cdf_type, "normal") == 0)
        norm = normal_gen_cdf(x, parameter_1, parameter_2);
    else if (strcmp(cdf_type, "uniform") == 0)
        norm = uniform_cdf(x, parameter_1, parameter_2);
    else if (strcmp(cdf_type, "tstudent") == 0)
        norm = tstudent_cdf(x, (int)parameter_1);
    else if (strcmp(cdf_type, "chi_squared") == 0)
        norm = chi_squared_cdf(x, (int)parameter_1);
    else if (strcmp(cdf_type, "exp") == 0)
        norm = exp_cdf(x, parameter_1);
    else if (strcmp(cdf_type, "cauchy") == 0)
        norm = Cauchy_cdf(x, parameter_1, parameter_2);
    else
    {
        fprintf(stderr, "Error! Allowed types of density are: ['normal', 'uniform', 'tstudent', 'chi_squared', 'exp', 'cauchy']!\n");
        exit(EXIT_FAILURE);
    }
    if (norm.args != NULL && norm.probs != NULL) {
        for (int i = 0; i < norm.size; i++) {

            PyList_Append(args_list, PyFloat_FromDouble(norm.args[i]));
            PyList_Append(probs_list, PyFloat_FromDouble(norm.probs[i]));
        }
    }

    PyDict_SetItemString(dict_res, "args", args_list);
    PyDict_SetItemString(dict_res, "probs", probs_list);
    PyDict_SetItemString(dict_res, "cdf", PyFloat_FromDouble(norm.cdf));

    //free(norm.args);
    //free(norm.probs);
    return dict_res;
}


static PyObject* percentile_cdf_py(PyObject* self, PyObject* args, PyObject* kwargs)
{
    double x;
    char* cdf_type = NULL;
    double parameter_1 = 0.0;
    double parameter_2 = 0.0;
    static char* keywords[] = { "x", "cdf_type", "parameter_1", "parameter_2", NULL };

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "ds|dd", &keywords, &x, &cdf_type, &parameter_1, &parameter_2))
        return NULL;

    if (strcmp(cdf_type, "normal") != 0 && strcmp(cdf_type, "cauchy") != 0 && strcmp(cdf_type, "uniform") != 0 &&
        strcmp(cdf_type, "tstudent") != 0 && strcmp(cdf_type, "chi_squared") != 0 && strcmp(cdf_type, "exp") != 0)
    {
        PyErr_SetString(PyExc_ValueError, "Unknown density type");
        return NULL;
    }

    if (strcmp(cdf_type, "normal") == 0)
        return Py_BuildValue("d", percentile_normal_gen_cdf(x, parameter_1, parameter_2));
    else if (strcmp(cdf_type, "uniform") == 0)
        return Py_BuildValue("d", percentile_uniform_cdf(x, parameter_1, parameter_2));
    else if (strcmp(cdf_type, "tstudent") == 0)
        return Py_BuildValue("d", percentile_tstudent_cdf(x, (int)parameter_1));
    else if (strcmp(cdf_type, "chi_squared") == 0)
        return Py_BuildValue("d", percentile_chi_squared_cdf(x, (int)parameter_1));
    else if (strcmp(cdf_type, "exp") == 0)
        return Py_BuildValue("d", percentile_exp_cdf(x, parameter_1));
    else if (strcmp(cdf_type, "cauchy") == 0)
        return Py_BuildValue("d", percentile_Cauchy_cdf(x, parameter_1, parameter_2));
    else
    {
        fprintf(stderr, "Error! Allowed types of density are: ['normal', 'uniform', 'tstudent', 'chi_squared', 'exp', 'cauchy']!\n");
        exit(EXIT_FAILURE);
    }
}


static PyObject* fib_py(PyObject* self, PyObject* args)
{
    PyObject* result = NULL;
    long n;
    if (!PyArg_ParseTuple(args, "l", &n))
        return NULL;
    result = Py_BuildValue("L", fib((unsigned int)n));

    return result;
}


static PyObject* normal_cdf_py(PyObject* self, PyObject* args)
{
    double x;
    if (!PyArg_ParseTuple(args, "d", &x))
        return NULL;

    struct cdf_data  norm = normal_cdf(x);

    PyObject* dict_res = PyDict_New();
    PyObject* args_list = PyList_New(0);
    PyObject* probs_list = PyList_New(0);

    if (norm.args != NULL && norm.probs != NULL) {
        for (int i = 0; i < norm.size; i++) {

            PyList_Append(args_list, PyFloat_FromDouble(norm.args[i]));
            PyList_Append(probs_list, PyFloat_FromDouble(norm.probs[i]));
        }
    }

    PyDict_SetItemString(dict_res, "args", args_list);
    PyDict_SetItemString(dict_res, "probs", probs_list);
    PyDict_SetItemString(dict_res, "cdf", PyFloat_FromDouble(norm.cdf));

    //free(norm.args);
    //free(norm.probs);
    return dict_res;
}


static PyObject* percentile_normal_cdf_py(PyObject* self, PyObject* args)
{
    double x;
    if (!PyArg_ParseTuple(args, "d", &x))
        return NULL;
    return Py_BuildValue("d", percentile_normal_cdf(x));

}



static char fib_doc[] = "Functon fib takes non-negative integer 'n'. Function fib returns n-th element of Fibonacci sequence.";
static char normal_cdf_doc[] = "Function takes floating, real argument x.\n Function returns arguments from -3 to x with step 0.01, correspondinf cdf values, and value of cdf corresponding presicely to x";
static char normal_density_fun_doc[] = "Function takes real argument x.\n Returns value of normal distribution N(0,1).";
static char percentile_normal_cdf_doc[] = "Takes float number from the interval (0,1).\n Returns percentile of the normal distribution.";
static char density_doc[] = "Function takes the following arguemnts:\n 'x' - floating number,\n 'density_type' - string, \n + parameters of choosen density (parameters have no default values set).\n Function returns value of chosen denisty evaluated at 'x' with provided parameters.";
static char cdf_doc[] = "Function take floating, real argument x, string name of cdf type and parameters of this function as floating numbers.\n  Function returns arguments with step 0.01, correspondinf cdf values, and value of cdf corresponding presicely to x";
static char percentile_cdf_doc[] = "Takes float number from the interval (0,1), string name of cdf type, and its parameters as floating numbers.\n Returns percentile of chosen distribution.";



static PyMethodDef math_functions_module_methods[] = {
    {"fibonacci", (PyCFunction)fib_py, METH_VARARGS, fib_doc},
    {"normal_cdf", (PyCFunction)normal_cdf_py, METH_VARARGS, normal_cdf_doc},
    {"dnormal", (PyCFunction)normal_density_py, METH_VARARGS, normal_density_fun_doc},
    {"percentile_normal", (PyCFunction)percentile_normal_cdf_py, METH_VARARGS, percentile_normal_cdf_doc},
    {"density", (PyCFunction)density_py, METH_VARARGS | METH_KEYWORDS, density_doc},
    {"cdf", (PyCFunction)cdf_py, METH_VARARGS | METH_KEYWORDS, cdf_doc},
    {"percentile_cdf", (PyCFunction)percentile_cdf_py, METH_VARARGS | METH_KEYWORDS, percentile_cdf_doc},
    {NULL, NULL, 0, NULL} };


static struct PyModuleDef math_functions_module_definition = {
    PyModuleDef_HEAD_INIT,
    "cmathfunctions",
    "C Extension of Python Code.",
    -1,
    math_functions_module_methods
};


PyMODINIT_FUNC PyInit_cmathfunctions(void) {
    return PyModule_Create(&math_functions_module_definition);
}


int main(void)
{
    return 0;
}
