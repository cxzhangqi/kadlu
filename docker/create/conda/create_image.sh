
# copy package, setup and environment file
cp -r ../../../kadlu/ .
cp ../../../setup.py .
cp ../../../environment.yml .

# build image
docker build --tag=kadlu_conda_test1 .

# tag image
docker tag kadlu_conda_test1 oliskir/kadlu:conda_test1

# push image to repository
docker push oliskir/kadlu:conda_test1

# clean
rm -rf kadlu
rm -rf setup.py
rm -rf environment.yml