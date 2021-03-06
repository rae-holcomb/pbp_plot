from setuptools import setup, find_packages

with open("README.md", "r") as readme_file:
    readme = readme_file.read()

requirements = ["lightkurve>=1"] # "ipython>=6", 

setup(
    name="pbp_plot",
    version="0.0.1",
    author="Rae Holcomb",
    author_email="raeholcomb19@gmail.com",
    description="Create pixel-by-pixel plots from TESS full frame images.",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/pbp_plot/homepage/",
    packages=find_packages(),
    install_requires=requirements,
    classifiers=[
        "Programming Language :: Python :: 3.6",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    ],
)
