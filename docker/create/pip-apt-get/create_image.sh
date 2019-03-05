
# copy package and requirements
cp -r ../../../kadlu/ .
cp ../../../setup.py .

# build image
docker build --tag=kadlu_test1 .

# tag image
docker tag kadlu_test1 oliskir/kadlu:test1

# push image to repository
docker push oliskir/kadlu:test1

# clean
rm -rf kadlu
rm -rf setup.py