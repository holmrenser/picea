Deploying to pypi:

```
poetry version <major,minor,patch>
poetry run coverage run
poetry run coverage report
poetry build
poetry deploy
```
