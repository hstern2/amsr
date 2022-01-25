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
# FromMol(m) takes an RDKit Mol m and returns an AMSR string
# ToSmiles(s) takes an AMSR string s and returns SMILES
# FromSmiles(s) takes a SMILES string s and returns AMSR
```

## Description

See [documentation](https://hstern2.github.io/amsr/) for
the structures produced by the examples below.

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
can still make bonds, according to its valence. Hydrogens may be added
explicitly like any other atom.  In the following
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


If no atoms are available to bond, a new molecule is formed.
```py
ToMol("CFFFFCO")
# carbon tetrafluoride and methanol
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
<path class='bond-0 atom-0 atom-1' d='M 100.0,79.3 L 127.4,79.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 127.4,79.3 L 154.8,79.3' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 100.0,79.3 L 72.6,79.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 72.6,79.3 L 45.2,79.3' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 100.0,79.3 L 100.0,53.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 100.0,53.7 L 100.0,28.0' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-0 atom-4' d='M 100.0,79.3 L 100.0,104.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-0 atom-4' d='M 100.0,104.9 L 100.0,130.6' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-5 atom-6' d='M 69.2,181.8 L 95.2,181.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-5 atom-6' d='M 95.2,181.8 L 121.2,181.8' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-1' d='M 156.3 70.6
L 166.7 70.6
L 166.7 72.6
L 158.7 72.6
L 158.7 77.9
L 165.8 77.9
L 165.8 79.9
L 158.7 79.9
L 158.7 88.0
L 156.3 88.0
L 156.3 70.6
' fill='#33CCCC'/>
<path  class='atom-2' d='M 33.3 70.6
L 43.7 70.6
L 43.7 72.6
L 35.6 72.6
L 35.6 77.9
L 42.8 77.9
L 42.8 79.9
L 35.6 79.9
L 35.6 88.0
L 33.3 88.0
L 33.3 70.6
' fill='#33CCCC'/>
<path  class='atom-3' d='M 94.8 9.1
L 105.2 9.1
L 105.2 11.1
L 97.2 11.1
L 97.2 16.4
L 104.3 16.4
L 104.3 18.4
L 97.2 18.4
L 97.2 26.5
L 94.8 26.5
L 94.8 9.1
' fill='#33CCCC'/>
<path  class='atom-4' d='M 94.8 132.1
L 105.2 132.1
L 105.2 134.1
L 97.2 134.1
L 97.2 139.4
L 104.3 139.4
L 104.3 141.4
L 97.2 141.4
L 97.2 149.5
L 94.8 149.5
L 94.8 132.1
' fill='#33CCCC'/>
<path  class='atom-6' d='M 122.8 181.9
Q 122.8 177.7, 124.8 175.4
Q 126.9 173.0, 130.8 173.0
Q 134.6 173.0, 136.7 175.4
Q 138.8 177.7, 138.8 181.9
Q 138.8 186.1, 136.7 188.5
Q 134.6 190.9, 130.8 190.9
Q 126.9 190.9, 124.8 188.5
Q 122.8 186.1, 122.8 181.9
M 130.8 188.9
Q 133.4 188.9, 134.8 187.2
Q 136.3 185.4, 136.3 181.9
Q 136.3 178.5, 134.8 176.7
Q 133.4 175.0, 130.8 175.0
Q 128.1 175.0, 126.6 176.7
Q 125.2 178.4, 125.2 181.9
Q 125.2 185.4, 126.6 187.2
Q 128.1 188.9, 130.8 188.9
' fill='#FF0000'/>
<path  class='atom-6' d='M 140.8 173.2
L 143.2 173.2
L 143.2 180.6
L 152.1 180.6
L 152.1 173.2
L 154.5 173.2
L 154.5 190.6
L 152.1 190.6
L 152.1 182.6
L 143.2 182.6
L 143.2 190.6
L 140.8 190.6
L 140.8 173.2
' fill='#FF0000'/>
</svg>

</div>



### Branches
Branches are formed automatically when atoms can no longer
make bonds.  They can also be made by "capping" or
"saturating" an atom with hydrogens, using a period.
The capping hydrogens are applied to the most
recently-added atom which can still make bonds.  In
the example below, the dot caps the third carbon,
so that the fourth is added to the second rather than the third,
to form isobutane rather than *n*-butane.
```py
ToMol("CCC.C")
# isobutane
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
<path class='bond-0 atom-0 atom-1' d='M 9.1,178.7 L 100.0,126.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 100.0,126.2 L 190.9,178.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-1 atom-3' d='M 100.0,126.2 L 100.0,21.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
</svg>

</div>



### Small rings
Small rings are denoted by one of the digits `3,4,5,6` indicating
the size of the ring.
A new bond is formed between atoms *i* and *j* where
*i* and *j* are the most recently-added atoms which
can still make bonds and when bonded will form a ring of that size.
```py
ToMol("CCO3")
# oxirane
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
<path class='bond-0 atom-0 atom-1' d='M 178.8,93.3 L 33.0,9.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-0' d='M 50.2,167.6 L 114.5,130.4' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-0' d='M 114.5,130.4 L 178.8,93.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 33.0,9.1 L 33.0,84.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 33.0,84.0 L 33.0,159.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-2' d='M 20.0 177.6
Q 20.0 170.8, 23.3 167.0
Q 26.7 163.2, 33.0 163.2
Q 39.3 163.2, 42.6 167.0
Q 46.0 170.8, 46.0 177.6
Q 46.0 184.5, 42.6 188.4
Q 39.2 192.3, 33.0 192.3
Q 26.7 192.3, 23.3 188.4
Q 20.0 184.5, 20.0 177.6
M 33.0 189.1
Q 37.3 189.1, 39.6 186.2
Q 42.0 183.3, 42.0 177.6
Q 42.0 172.0, 39.6 169.2
Q 37.3 166.4, 33.0 166.4
Q 28.7 166.4, 26.3 169.2
Q 24.0 172.0, 24.0 177.6
Q 24.0 183.3, 26.3 186.2
Q 28.7 189.1, 33.0 189.1
' fill='#FF0000'/>
</svg>

</div>


```py
ToMol("CCCCCC6")
# cyclohexane
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
<path class='bond-0 atom-0 atom-1' d='M 190.9,100.0 L 145.5,178.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-0' d='M 145.5,21.3 L 190.9,100.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 145.5,178.7 L 54.5,178.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 54.5,178.7 L 9.1,100.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 9.1,100.0 L 54.5,21.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 54.5,21.3 L 145.5,21.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
</svg>

</div>



### Large rings
Rings with seven or more members are denoted by a
sequence of digits `3,4,5,6` which are added to give the size of the ring.
```py
ToMol("CCCCCCC43")
# cycloheptane
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
<path class='bond-0 atom-0 atom-1' d='M 188.6,100.0 L 153.5,27.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-6 atom-0' d='M 153.5,172.9 L 188.6,100.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 153.5,27.1 L 74.6,9.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 74.6,9.1 L 11.4,59.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 11.4,59.5 L 11.4,140.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 11.4,140.5 L 74.6,190.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-6' d='M 74.6,190.9 L 153.5,172.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
</svg>

</div>


```py
ToMol("CCCCCCCCCC55")
# cyclodecane
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
<path class='bond-0 atom-0 atom-1' d='M 190.9,100.0 L 173.5,153.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-9 atom-0' d='M 173.5,46.6 L 190.9,100.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 173.5,153.4 L 128.1,186.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 128.1,186.5 L 71.9,186.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 71.9,186.5 L 26.5,153.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 26.5,153.4 L 9.1,100.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-6' d='M 9.1,100.0 L 26.5,46.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-6 atom-7' d='M 26.5,46.6 L 71.9,13.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-7 atom-8' d='M 71.9,13.5 L 128.1,13.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-8 atom-9' d='M 128.1,13.5 L 173.5,46.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
</svg>

</div>


If any characters appear between digits, multiple rings are formed
rather than a single large ring.
Whitespace may be used if this is desired (without adding any additional atoms).
```py
ToMol("CCCCCCCCCC5 5")
# 1-ethyl-1,2,3,3a,4,5,6,6a-octahydropentalene
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
<path class='bond-0 atom-0 atom-1' d='M 190.9,80.1 L 150.7,65.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 150.7,65.2 L 117.8,92.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 117.8,92.5 L 120.5,135.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-2' d='M 76.2,81.9 L 117.8,92.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 120.5,135.3 L 80.7,151.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 80.7,151.1 L 53.3,118.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-6' d='M 53.3,118.1 L 11.8,107.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-9 atom-5' d='M 76.2,81.9 L 53.3,118.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-6 atom-7' d='M 11.8,107.5 L 9.1,64.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-7 atom-8' d='M 9.1,64.7 L 48.9,48.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-8 atom-9' d='M 48.9,48.9 L 76.2,81.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
</svg>

</div>



## Developing
This repo uses pre-commit, so after cloning run `pip install -r requirements.txt` and `pre-commit install` prior to committing.
