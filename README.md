# selenzy-wrapper - rpSBML wrapper for selenzy tool

<!-- ## Install
### From Conda
```sh
[sudo] conda install -c brsynth -c chemlite
```

## Use
### Compound
```python
from chemlite import Compound

c = Compound(id='test_cmpd')
```
The code above creates an empty compound. The following fields can be filled and accessed either at build time or later on:
- smiles
- inchi
- inchikey
- formula
- name
- infos

### Reaction
```python
from chemlite import Reaction

r = Reaction(id='test_rxn')
```
The code above creates an empty reaction. The following fields can be filled and accessed either at build time or later on:
- ec_numbers
- reactants
- products
- infos

The following methods are also available:
- `get_smiles()`
- `add_reactant()`
- `add_product()`


### Pathway
```python
from chemlite import Pathway

p = Pathway(id='test_path')
```
The code above creates an empty reaction. The following fields can be filled and accessed either at build time or later on:
- id
- species
- reactions

The following methods are also available:
- `add_compound()`
- `add_reaction()`
- `del_reaction()`
- `Pathway.net_reaction()`


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
