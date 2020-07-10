Deploying to pypi:
```
rm -r dist
python setup.py sdist bdist_wheel
twine upload dist/*
```