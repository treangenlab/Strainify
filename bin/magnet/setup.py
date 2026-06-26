from setuptools import setup

setup(
    name="magnet",
    version="0.1.0",
    py_modules=["magnet"],
    packages=["utils"],
    install_requires=[
        "pandas",
        "biopython",
        "ete3",
        "scikit-learn",
        "numpy",
    ],
    entry_points={
        "console_scripts": [
            # calls main() in magnet.py
            "magnet=magnet:main",
        ],
    },
)

