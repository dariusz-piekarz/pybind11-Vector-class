from setuptools import setup, Extension


sfc_module = Extension('cmathfunctions', sources = ['percentile_py_d.c'])

setup(
    name = 'cmathfunctions',
    version = '1.0',
    description = 'Python functions in C',
    ext_modules = [sfc_module],
    compiler_directives=dict(always_allow_keywords=True)
)
