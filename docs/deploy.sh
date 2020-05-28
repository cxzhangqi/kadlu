#!/bin/bash

# Note: the host kadlu_docs must be configured in ~/.ssh/config 

# kadlu source directory
if [ ! -d $HOME/kadlu ]; then 
   echo 'Enter path to kadlu folder, e.g., /home/user/'
   read path
   KADLUDIR="$path/kadlu"
else
   KADLUDIR="$HOME/kadlu"
fi

# export the version as an environment variable
# this will be used in setup.py and documentation conf.py
echo 'Enter the version number'
read KADLUVERSION
export KADLUVERSION=$KADLUVERSION

# edit the changelog
vim $KADLUDIR/docs/source/versions/changelog.rst

# update sphinx, install sphinx theme, build the docs
python3 -m pip install --upgrade sphinx sphinx_rtd_theme $KADLUDIR/docs/sphinx_mer_rtd_theme-0.4.3.dev0.tar.gz
make -C $KADLUDIR/docs html

# copy the html content to the docs server
ssh kadlu_docs << EOF
sudo mkdir /var/www/html/kadlu2
sudo chown -h -R ubuntu:ubuntu /var/www/html/kadlu2
sudo chown -h -R ubuntu:ubuntu /var/www/html/kadlu
EOF
scp -r $KADLUDIR/docs/build/html/* kadlu_docs:/var/www/html/kadlu2
ssh kadlu_docs << EOF
sudo rm -rf /var/www/html/kadlu
sudo mv /var/www/html/kadlu2 /var/www/html/kadlu
sudo chown -h -R www-data:www-data /var/www/html/kadlu
sudo chmod -R g+w /var/www/html/kadlu
EOF

# build the package and upload to pypi
cd $KADLUDIR
python3 $KADLUDIR/setup.py sdist bdist_wheel
python3 -m twine upload --repository pypi $KADLUDIR/dist/*

