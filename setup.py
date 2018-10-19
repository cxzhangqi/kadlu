from setuptools import setup

setup(name='mafmlib',
      version='0.0.1',
      description="Python toolbox for the MERIDIAN Acoustic Field Model",
      url='https://data.meridian.cs.dal.ca/gitlab/data_analytics_dal/packages/mafmlib',
      author='Mark Thomas, Oliver Kirsebom',
      author_email='markthomas@dal.ca, oliver.kirsebom@dal.ca',
      license='GNU General Public License v3.0',
      packages=['mafmlib'],
      install_requires=[
          'numpy',
          'arlpy',
          ],
      setup_requires=['pytest-runner', ],
      tests_require=['pytest', ],
      include_package_data=True,
      zip_safe=False)
