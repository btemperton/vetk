from setuptools import setup


def readme():
    with open('README.rst') as f:
        return f.read()


setup(name='ViralEcologyToolkit',
      version='0.1',
      description='Collection of tools for the Temperton lab',
      url='https://git.exeter.ac.uk/bt273/ViralEcologyToolkit',
      author='Ben Temperton',
      author_email='b.temperton@exeter.ac.uk',
      license='MIT',
      packages=['ViralEcologyToolkit'],
      entry_points={
          'console_scripts': ['create_reticulate_network=ReticulateNetwork.command_line:main'],
      },
      install_requires=['numpy', 'pandas', 'scipy', 'biopython', 'rpy2'],
      zip_safe=False)
