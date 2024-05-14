

## Dependencies 
We use [`poetry`](https://python-poetry.org) to version dependencies and publish to [PyPi](https://pypi.org).

Once you have `poetry`, you may initialize a new project by entering its root folder and running:
```bash
poetry init
```
Then you may add dependencies simple by running:
```bash
poetry add numpy
```
where `numpy` can be any dependency required. 

Then, you may install the dependencies by:
```bash
poetry install
```
and build it :

```bash
poetry build
```

Finally, once a new version is ready, one may release it on PyPi by running:
```bash
poetry publish
```
Yet, credentials are required for the last step.


### Notes on dependencies

We need `seaborn` at least `0.13.0` otherwise annotation in heatmaps fails.

```bash
python -m pip install seaborn==0.13.0
```

In case we need a `requirements.txt` file, we only have to run:
```bash
poetry export --without-hashes --format=requirements.txt > requirements.txt
```


