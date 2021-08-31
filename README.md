# SimBa2

SimBa2 is a Python package for predicting changes in protein stability upon mutation using the prediction methods SimBa-IB and SimBa-SYM. The program takes a PDB file as input, and outputs the predicted stability changes (in kcal/mol) for all possible single amino acid changes.

Negative ddG values signify destabilizing mutations.

The program can be run from the command line or used as a Python module.

## Installation

```
python -m pip install https://github.com/ktbaek/SimBa2/archive/main.tar.gz
```

## Usage

### Run from command line

The program takes a single PDB code or a local PDB file as input, and outputs two .csv files (SimBa-IB_<pdb-code>.csv and SimBa-SYM_<pdb-code>.csv) with predicted ddG values in kcal/mol for a saturated mutagenesis of the protein.

#### Examples

Analyze PDB structure 3BDN:
```
python -m simba2 3bdn
```
Keep downloaded PDB file in the current directory:
```
python -m simba2 3bdn --keep
```
Keep downloaded PDB file in another directory:
```
python -m simba2 3bdn --dir ../my_directory --keep
```
Analyze local PDB file of 3BDN:
```
python -m simba2 3bdn --file 3bdn.pdb
```
Help:
```
python -m simba2 --help
```

### Used as Python module

The package can also be used as a module in a Python script to analyze local PDB files.

Example:

```
from simba2 import methods as sb
sb.simba2_predict('3BDN', '3bdn.pdb')
```

### Protein structure file format

The package so far only works with structures that are available in PDB format, not mmCIF format.

## Output

The resulting .csv files (or dataframe when used in a Python script) have mostly self-explanatory column names.

The column named 'mean_RSA' shows the mean RSA value across chains if the structure is a homooligomer.

The column named 'Gene' shows the chain-assigned gene name as it appears in the header section of the PDB file. This can be useful when analyzing heterooligomeric proteins. Sometimes this column is empty.

## Citing SimBa2
Please acknowledge the use of SimBa2 by citing the following publication:

Bæk, K. T. and K. P. Kepp (2021) Dataset and fitting dependencies when estimating protein mutant stability: Towards simple, balanced, interpretable models

Journal: TBA

DOI: TBA

## Copyright and license
The MIT License (MIT)

Copyright (c) 2021 Kristoffer Torbjørn Bæk, Kasper Planeta Kepp, Technical University of Denmark

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
