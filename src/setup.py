from distutils.core import setup, Extension

cpp_args = ['-std=c++17', '-lpthread']

ext_modules = [
    Extension(
        'prime_gen',
        ['prime_gen.cpp', "big_integer.cpp"],
        include_dirs=['pybind11/include'],
        language='c++',
        extra_compile_args = cpp_args,
    ),
]

setup(
    name='prime_gen',
    version='1.0',
    author='Dan Strano',
    author_email='dan@unitary.fund',
    description='pybind11 plugin to generate prime numbers',
    ext_modules=ext_modules,
)
