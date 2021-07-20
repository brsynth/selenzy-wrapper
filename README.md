# selenzy-wrapper - rpSBML wrapper for selenzy tool

## Install
### From Conda
```sh
[sudo] conda install -c brsynth -c conda-forge -c bioconda selenzy_wrapper

```

## Use
```python
from selenzy_wrapper import selenzy_pathway

pathway = rpPathway.from_rpSBML(infile='tests/data/pathway.xml')

selenzy_pathway(pathway=pathway)

pathway.to_rpSBML().write_to_file('selenzy.xml')
```

## Tests
Please follow instructions below ti run tests:
```
cd tests
pytest -v
```
For further tests and development tools, a CI toolkit is provided in `ci` folder (see [ci/README.md](ci/README.md)). -->


## Authors

* **Joan Hérisson**

## Acknowledgments

* Thomas Duigou


## Licence
selenzy-wrapper is released under the MIT licence. See the LICENCE file for details.
