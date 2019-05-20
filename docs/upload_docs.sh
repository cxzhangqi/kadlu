#!/bin/bash

# prompt user for username to MERIDIAN's docs server
echo 'Enter username'
read user

# copy html folder
scp -P 10022 -r build/html $user@206.12.88.81:/var/www/html/kadlu2

# replace old folder on server, and set permissions
ssh -p 10022 $user@206.12.88.81 << EOF
cd /var/www/html/
rm -rf kadlu
mv kadlu2 kadlu
sudo chown -h -R www-data:www-data kadlu
sudo chmod -R g+w kadlu
EOF

