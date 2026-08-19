import sys

from setuptools import Extension, setup

extra_link_args = ["-framework", "Accelerate"] if sys.platform == "darwin" else []
libraries = [] if sys.platform == "darwin" else ["m"]

setup(
    ext_modules=[
        Extension(
            "amsr.src.conf_util",
            sources=["amsr/src/conf_util.c", "amsr/src/lbfgs.c"],
            include_dirs=["amsr/src"],
            extra_compile_args=["-O3"],
            extra_link_args=extra_link_args,
            libraries=libraries,
        )
    ]
)
