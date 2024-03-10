from setuptools import setup, Extension
import pybind11

cpp_args = ['-std=c++11', '-stdlib=libc++', '-mmacosx-version-min=10.7']

sfc_module = Extension(
    'math_module',
    sources=['vector_d.cpp'],
    include_dirs=[pybind11.get_include()],
    language='c++',
    extra_compile_args=cpp_args,
    )

setup(
    name='math_module',
    version='1.0',
    description='Fibonacci and Vector C++ class.',
    ext_modules=[sfc_module],
)


