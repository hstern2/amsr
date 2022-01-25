# AMSR
**A**nother **M**olecular **S**tring **R**epresentation,
inspired by [H. Hiz "A Linearization of Chemical Graphs,"](https://pubs.acs.org/doi/10.1021/c160014a015) [SMILES](https://pubs.acs.org/doi/10.1021/ci00057a005), [DeepSMILES](https://github.com/baoilleach/deepsmiles), and [SELFIES](https://github.com/aspuru-guzik-group/selfies).

## Installing
```
pip install git+git://github.com/hstern2/amsr.git
```

## Usage
```py
from amsr import ToMol, FromMol, ToSmiles, FromSmiles
# ToMol(s) takes an AMSR string s and returns an RDKit Mol
# FromMol(m) takes an RDKit Mol m and returns a string
# ToSmiles(s) takes an AMSR string s and returns SMILES
# FromSmiles(s) takes a SMILES string s and returns AMSR
```

## Description

### Atoms
As in SMILES, atoms are represented by their symbol
enclosed in square brackets. If the atom is designated by a one-letter symbol,
brackets may be omitted.  All atoms are assumed to have a fixed
valence.  If an atom makes fewer bonds than its valence, hydrogens are assumed.

```py
ToMol("C")
# methane
```

<div>
<?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='200px' height='200px' viewBox='0 0 200 200'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='200.0' height='200.0' x='0.0' y='0.0'> </rect>
<path  class='atom-0' d='M 89.1 100.1
Q 89.1 93.1, 92.4 89.4
Q 95.7 85.7, 102.0 85.7
Q 107.8 85.7, 110.9 89.8
L 108.3 92.0
Q 106.0 89.0, 102.0 89.0
Q 97.7 89.0, 95.4 91.8
Q 93.2 94.7, 93.2 100.1
Q 93.2 105.7, 95.5 108.6
Q 97.8 111.5, 102.4 111.5
Q 105.5 111.5, 109.2 109.6
L 110.3 112.6
Q 108.8 113.6, 106.6 114.1
Q 104.3 114.7, 101.8 114.7
Q 95.7 114.7, 92.4 110.9
Q 89.1 107.2, 89.1 100.1
' fill='#000000'/>
<path  class='atom-0' d='M 114.3 86.0
L 118.2 86.0
L 118.2 98.0
L 132.6 98.0
L 132.6 86.0
L 136.5 86.0
L 136.5 114.3
L 132.6 114.3
L 132.6 101.2
L 118.2 101.2
L 118.2 114.3
L 114.3 114.3
L 114.3 86.0
' fill='#000000'/>
<path  class='atom-0' d='M 149.4 122.1
L 151.7 122.1
L 151.7 124.2
L 149.4 124.2
L 149.4 128.5
L 146.9 128.5
L 146.9 124.2
L 137.3 124.2
L 137.3 122.5
L 145.5 109.8
L 149.4 109.8
L 149.4 122.1
M 140.3 122.1
L 146.9 122.1
L 146.9 111.5
L 140.3 122.1
' fill='#000000'/>
</svg>

</div>


```py
ToMol("O")
# water
```

<div>
<?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='200px' height='200px' viewBox='0 0 200 200'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='200.0' height='200.0' x='0.0' y='0.0'> </rect>
<path  class='atom-0' d='M 48.9 86.0
L 52.8 86.0
L 52.8 98.0
L 67.2 98.0
L 67.2 86.0
L 71.1 86.0
L 71.1 114.3
L 67.2 114.3
L 67.2 101.2
L 52.8 101.2
L 52.8 114.3
L 48.9 114.3
L 48.9 86.0
' fill='#FF0000'/>
<path  class='atom-0' d='M 72.5 113.3
Q 73.1 111.6, 74.8 110.6
Q 76.4 109.6, 78.7 109.6
Q 81.5 109.6, 83.1 111.1
Q 84.7 112.6, 84.7 115.4
Q 84.7 118.1, 82.6 120.7
Q 80.6 123.3, 76.4 126.4
L 85.0 126.4
L 85.0 128.5
L 72.4 128.5
L 72.4 126.7
Q 75.9 124.2, 78.0 122.4
Q 80.0 120.5, 81.0 118.9
Q 82.0 117.2, 82.0 115.5
Q 82.0 113.7, 81.1 112.7
Q 80.2 111.7, 78.7 111.7
Q 77.2 111.7, 76.2 112.3
Q 75.2 112.9, 74.5 114.3
L 72.5 113.3
' fill='#FF0000'/>
<path  class='atom-0' d='M 87.0 100.1
Q 87.0 93.3, 90.4 89.5
Q 93.7 85.7, 100.0 85.7
Q 106.3 85.7, 109.6 89.5
Q 113.0 93.3, 113.0 100.1
Q 113.0 107.0, 109.6 110.9
Q 106.2 114.8, 100.0 114.8
Q 93.8 114.8, 90.4 110.9
Q 87.0 107.0, 87.0 100.1
M 100.0 111.6
Q 104.3 111.6, 106.6 108.7
Q 109.0 105.8, 109.0 100.1
Q 109.0 94.5, 106.6 91.7
Q 104.3 88.9, 100.0 88.9
Q 95.7 88.9, 93.3 91.7
Q 91.0 94.5, 91.0 100.1
Q 91.0 105.8, 93.3 108.7
Q 95.7 111.6, 100.0 111.6
' fill='#FF0000'/>
</svg>

</div>


```py
ToMol("[Cl]")
# hydrochloric acid
```

<div>
<?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='200px' height='200px' viewBox='0 0 200 200'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='200.0' height='200.0' x='0.0' y='0.0'> </rect>
<path  class='atom-0' d='M 48.7 86.9
L 52.6 86.9
L 52.6 98.9
L 67.0 98.9
L 67.0 86.9
L 70.9 86.9
L 70.9 115.2
L 67.0 115.2
L 67.0 102.1
L 52.6 102.1
L 52.6 115.2
L 48.7 115.2
L 48.7 86.9
' fill='#00CC00'/>
<path  class='atom-0' d='M 72.9 101.0
Q 72.9 93.9, 76.2 90.3
Q 79.5 86.5, 85.8 86.5
Q 91.6 86.5, 94.7 90.7
L 92.1 92.8
Q 89.8 89.8, 85.8 89.8
Q 81.5 89.8, 79.2 92.7
Q 77.0 95.5, 77.0 101.0
Q 77.0 106.6, 79.3 109.5
Q 81.6 112.3, 86.2 112.3
Q 89.3 112.3, 93.0 110.5
L 94.1 113.5
Q 92.6 114.4, 90.4 115.0
Q 88.1 115.5, 85.6 115.5
Q 79.5 115.5, 76.2 111.8
Q 72.9 108.0, 72.9 101.0
' fill='#00CC00'/>
<path  class='atom-0' d='M 98.2 84.8
L 101.8 84.8
L 101.8 115.2
L 98.2 115.2
L 98.2 84.8
' fill='#00CC00'/>
</svg>

</div>



### Chains
Each atom in a chain is bonded to the most recently added atom that
can still make bonds, according to its valence. In the following
example, the fluorines are added to the second carbon.  The chlorine
is then added to the first carbon, since the second can no longer bond.
```py
ToMol("CCFFF[Cl]")
# 2-chloro-1,1,1-trifluoroethane
```

<div>
<?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='200px' height='200px' viewBox='0 0 200 200'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='200.0' height='200.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 79.1,130.8 L 132.4,100.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-0 atom-5' d='M 79.1,130.8 L 53.8,116.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-0 atom-5' d='M 53.8,116.2 L 28.5,101.5' style='fill:none;fill-rule:evenodd;stroke:#00CC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 132.4,100.0 L 155.7,86.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 155.7,86.6 L 179.0,73.1' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-1 atom-3' d='M 132.4,100.0 L 144.8,121.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-1 atom-3' d='M 144.8,121.5 L 157.3,143.1' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-1 atom-4' d='M 132.4,100.0 L 120.0,78.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-1 atom-4' d='M 120.0,78.5 L 107.6,56.9' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-2' d='M 180.5 60.5
L 190.9 60.5
L 190.9 62.5
L 182.9 62.5
L 182.9 67.8
L 190.0 67.8
L 190.0 69.8
L 182.9 69.8
L 182.9 77.9
L 180.5 77.9
L 180.5 60.5
' fill='#33CCCC'/>
<path  class='atom-3' d='M 158.0 144.6
L 168.4 144.6
L 168.4 146.6
L 160.3 146.6
L 160.3 151.9
L 167.5 151.9
L 167.5 153.9
L 160.3 153.9
L 160.3 162.0
L 158.0 162.0
L 158.0 144.6
' fill='#33CCCC'/>
<path  class='atom-4' d='M 96.4 38.0
L 106.8 38.0
L 106.8 40.0
L 98.8 40.0
L 98.8 45.3
L 105.9 45.3
L 105.9 47.3
L 98.8 47.3
L 98.8 55.4
L 96.4 55.4
L 96.4 38.0
' fill='#33CCCC'/>
<path  class='atom-5' d='M 9.1 100.6
Q 9.1 96.3, 11.1 94.0
Q 13.2 91.7, 17.0 91.7
Q 20.6 91.7, 22.5 94.3
L 20.9 95.6
Q 19.5 93.7, 17.0 93.7
Q 14.4 93.7, 13.0 95.5
Q 11.6 97.3, 11.6 100.6
Q 11.6 104.1, 13.0 105.8
Q 14.5 107.6, 17.3 107.6
Q 19.2 107.6, 21.5 106.4
L 22.1 108.3
Q 21.2 108.9, 19.9 109.2
Q 18.5 109.6, 16.9 109.6
Q 13.2 109.6, 11.1 107.3
Q 9.1 104.9, 9.1 100.6
' fill='#00CC00'/>
<path  class='atom-5' d='M 24.7 90.7
L 26.9 90.7
L 26.9 109.3
L 24.7 109.3
L 24.7 90.7
' fill='#00CC00'/>
</svg>

</div>



## Developing
This repo uses pre-commit, so after cloning run `pip install -r requirements.txt` and `pre-commit install` prior to committing.
