from setuptools import setup, find_packages

setup(name="synthetic",
      packages=find_packages(),
      description="synthetic cluster injection to full line of sight image simulations",
      install_requires=['numpy', 'scipy', 'pandas', ],
      author="Tamas Norbert Varga",
      author_email="T.Varga@physik.lmu.de",
      version="0.1")