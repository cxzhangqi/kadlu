from setuptools import setup, find_packages

#
# create distribution and upload to pypi.org with:
#  
#   $ python setup.py sdist bdist_wheel
#   $ twine upload dist/*
#

setup(name='kadlu',
      version='0.0.1',
      description="MERIDIAN's Ocean Soundscape Toolbox (OST) written in Python",
      url='https://data.meridian.cs.dal.ca/gitlab/data_analytics_dal/packages/kadlu',
      author='Mark Thomas, Oliver Kirsebom',
      author_email='markthomas@dal.ca, oliver.kirsebom@dal.ca',
      license='GNU General Public License v3.0',
      packages=find_packages(),
      install_requires=[
          'numpy',
          'arlpy',
          ],
      setup_requires=['pytest-runner', ],
      tests_require=['pytest', ],
      include_package_data=True,
      zip_safe=False)
