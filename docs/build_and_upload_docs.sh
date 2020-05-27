#!/bin/bash

# Note: the host kadlu_docs must be configured in your ~/.ssh/config 

echo 'Enter path to kadlu folder, e.g., /home/user/'
read path

# kadlu source directory
KADLUDIR="$path/kadlu"

# update sphinx, install sphinx theme, build the docs
python3 -m pip install --upgrade pip
python3 -m pip install sphinx sphinx_rtd_theme
cd $KADLUDIR/docs 
python3 -m pip install sphinx_mer_rtd_theme-0.4.3.dev0.tar.gz
make html

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
