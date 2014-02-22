cd docs/doxygen
doxgyen Doxyfile
cd ..
make html
cd ..
cp -rf docs/build/html/. .
git add .
git add -f _images/* .
git add -f _images/math/* .
git add -f doxygen/html/* .
