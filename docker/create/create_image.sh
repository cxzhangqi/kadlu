
# copy package and requirements
cp -r ../../pyost/ .
cp ../../requirements.txt .
cp ../../setup.py .

# build image
docker build --tag=pyost_test1 .

# tag image
docker tag pyost_test1 oliskir/pyost:test1

# push image to repository
docker push oliskir/pyost:test1

# clean
rm -rf pyost
rm -rf requirements
rm -rf setup.py