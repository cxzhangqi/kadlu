from setuptools import setup, find_packages

# create distribution and upload to pypi.org with:
#   $ python setup.py sdist bdist_wheel
#   $ twine upload dist/*

setup(name='kadlu',
      version='1.0.0',
      description="MERIDIAN Python package for ocean ambient noise modelling",
      url='https://gitlab.meridian.cs.dal.ca/public_projects/kadlu',
      author='Matthew Smith, Oliver Kirsebom',
      author_email='matthew.smith@dal.ca, oliver.kirsebom@dal.ca',
      license='GNU General Public License v3.0',
      packages=find_packages(),
      install_requires=[
          'numpy',
          ],
      setup_requires=['pytest-runner', ],
      tests_require=['pytest', ],
      include_package_data=True,
      zip_safe=False)
