
# copy package, setup and environment file
cp -r ../../../kadlu/ .
cp ../../../setup.py .
cp ../../../environment.yml .
cp ../../../../meridian-rtd-theme/dist/sphinx_mer_rtd_theme-0.4.3.dev0.tar.gz .

# build image
docker build --tag=kadlu_v0.0.3 .

# tag image
docker tag kadlu_v0.0.3 oliskir/kadlu:v0.0.3

# push image to repository
docker push oliskir/kadlu:v0.0.3

# clean
rm -rf kadlu
rm -rf setup.py
rm -rf environment.yml
rm -rf sphinx_mer_rtd_theme-0.4.3.dev0.tar.gz