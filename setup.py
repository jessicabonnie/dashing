from setuptools import setup, Extension
import pybind11

ext_modules = [
    Extension(
        "dashing",
        ["src/bindings/dashing_bind.cpp"],
        include_dirs=[
            pybind11.get_include(),
            "src"
        ],
        language='c++'
    ),
]

setup(
    name="dashing",
    ext_modules=ext_modules,
    install_requires=['pybind11>=2.6.0'],
    python_requires=">=3.6"
)