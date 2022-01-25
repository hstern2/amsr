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

See [documentation](https://hstern2.github.io/amsr/)
for structures accompanying the following examples.

### Atoms
As in SMILES, atoms are represented by their symbol
enclosed in square brackets. If the atom is designated by a one-letter symbol,
brackets may be omitted.  All atoms are assumed to have a fixed
valence.  If an atom makes fewer bonds than its valence, hydrogens are assumed.



| molecule | AMSR |
| ---  | ---  |
| [methane](https://en.wikipedia.org/wiki/Methane) | C |

<div><?xml version='1.0' encoding='iso-8859-1'?>
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




| molecule | AMSR |
| ---  | ---  |
| [water](https://en.wikipedia.org/wiki/Water) | O |

<div><?xml version='1.0' encoding='iso-8859-1'?>
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




| molecule | AMSR |
| ---  | ---  |
| [hydrochloric acid](https://en.wikipedia.org/wiki/Hydrochloric_acid) | [Cl] |

<div><?xml version='1.0' encoding='iso-8859-1'?>
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


| molecule | AMSR |
| ---  | ---  |
| [2-chloro-1,1,1-trifluoroethane](https://pubchem.ncbi.nlm.nih.gov/compound/2-Chloro-1_1_1-trifluoroethane) | CCFFF[Cl] |

<div><?xml version='1.0' encoding='iso-8859-1'?>
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


| molecule | AMSR |
| ---  | ---  |
| [carbon tetrafluoride](https://en.wikipedia.org/wiki/Carbon_tetrafluoride) and [methanol](https://en.wikipedia.org/wiki/Methanol) | CFFFFCO |

<div><?xml version='1.0' encoding='iso-8859-1'?>
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
"saturating" an atom with hydrogens, using a period (`.`).
The capping hydrogens are applied to the most
recently-added atom which can still make bonds.  In
the example below, the dot caps the third carbon,
so that the fourth is added to the second rather than the third,
to form isobutane rather than *n*-butane.


| molecule | AMSR |
| ---  | ---  |
| [isobutane](https://en.wikipedia.org/wiki/Isobutane) | CCC.C |

<div><?xml version='1.0' encoding='iso-8859-1'?>
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


| molecule | AMSR |
| ---  | ---  |
| [oxirane](https://en.wikipedia.org/wiki/Ethylene_oxide) | CCO3 |

<div><?xml version='1.0' encoding='iso-8859-1'?>
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




| molecule | AMSR |
| ---  | ---  |
| [cyclohexane](https://en.wikipedia.org/wiki/Cyclohexane) | CCCCCC6 |

<div><?xml version='1.0' encoding='iso-8859-1'?>
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


| molecule | AMSR |
| ---  | ---  |
| [cycloheptane](https://en.wikipedia.org/wiki/Cycloheptane) | CCCCCCC43 |

<div><?xml version='1.0' encoding='iso-8859-1'?>
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




| molecule | AMSR |
| ---  | ---  |
| [cyclododecane](https://en.wikipedia.org/wiki/Cyclododecane) | CCCCCCCCCCCC66 |

<div><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='200px' height='200px' viewBox='0 0 200 200'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='200.0' height='200.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 190.9,100.0 L 178.7,145.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-11 atom-0' d='M 178.7,54.5 L 190.9,100.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 178.7,145.5 L 145.5,178.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 145.5,178.7 L 100.0,190.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 100.0,190.9 L 54.5,178.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 54.5,178.7 L 21.3,145.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-6' d='M 21.3,145.5 L 9.1,100.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-6 atom-7' d='M 9.1,100.0 L 21.3,54.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-7 atom-8' d='M 21.3,54.5 L 54.5,21.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-8 atom-9' d='M 54.5,21.3 L 100.0,9.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-9 atom-10' d='M 100.0,9.1 L 145.5,21.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-10 atom-11' d='M 145.5,21.3 L 178.7,54.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
</svg>
</div>



If any characters appear between digits, multiple rings are formed
rather than a single large ring.
Whitespace may be used if this is desired (without adding any additional atoms).


| molecule | AMSR |
| ---  | ---  |
| [1-ethyldecalin](https://pubchem.ncbi.nlm.nih.gov/compound/33053) | CCCCCCCCCCCC6 6 |

<div><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='200px' height='200px' viewBox='0 0 200 200'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='200.0' height='200.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 9.1,55.8 L 49.0,53.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 49.0,53.0 L 71.3,86.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 71.3,86.2 L 53.8,122.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-12 atom-11 atom-2' d='M 111.2,83.4 L 71.3,86.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 53.8,122.1 L 76.1,155.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 76.1,155.3 L 116.0,152.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-6' d='M 116.0,152.5 L 133.5,116.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-6 atom-7' d='M 133.5,116.6 L 173.4,113.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-11 atom-6' d='M 111.2,83.4 L 133.5,116.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-7 atom-8' d='M 173.4,113.8 L 190.9,77.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-8 atom-9' d='M 190.9,77.9 L 168.6,44.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-9 atom-10' d='M 168.6,44.7 L 128.7,47.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-10 atom-11' d='M 128.7,47.5 L 111.2,83.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
</svg>
</div>



### Double bonds (sp<sup>2</sup> centers)
Atoms making a single double bond are indicated by changing
the symbol to lowercase (note that lowercase
does not mean "aromatic"; merely, "atom having one
fewer neighbor than its valence.")  Double bonds are
assigned by a matching algorithm.  If a perfect
matching cannot be found (for instance, in the case
of an odd number of contiguous lowercase atoms) a
maximal matching is chosen, non-matched atoms
remain singly bonded, and hydrogens are added.


| molecule | AMSR |
| ---  | ---  |
| [formaldehyde](https://en.wikipedia.org/wiki/Formaldehyde) | co |

<div><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='200px' height='200px' viewBox='0 0 200 200'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='200.0' height='200.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 25.9,113.6 L 85.9,113.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 85.9,113.6 L 145.9,113.6' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 25.9,86.4 L 85.9,86.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 85.9,86.4 L 145.9,86.4' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-1' d='M 149.3 100.1
Q 149.3 93.3, 152.6 89.5
Q 156.0 85.7, 162.3 85.7
Q 168.6 85.7, 171.9 89.5
Q 175.3 93.3, 175.3 100.1
Q 175.3 107.0, 171.9 110.9
Q 168.5 114.8, 162.3 114.8
Q 156.0 114.8, 152.6 110.9
Q 149.3 107.0, 149.3 100.1
M 162.3 111.6
Q 166.6 111.6, 168.9 108.7
Q 171.3 105.8, 171.3 100.1
Q 171.3 94.5, 168.9 91.7
Q 166.6 88.9, 162.3 88.9
Q 158.0 88.9, 155.6 91.7
Q 153.3 94.5, 153.3 100.1
Q 153.3 105.8, 155.6 108.7
Q 158.0 111.6, 162.3 111.6
' fill='#FF0000'/>
</svg>
</div>




| molecule | AMSR |
| ---  | ---  |
| [benzene](https://en.wikipedia.org/wiki/Benzene) | cccccc6 |

<div><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='200px' height='200px' viewBox='0 0 200 200'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='200.0' height='200.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 190.9,100.0 L 145.5,178.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 168.3,102.7 L 136.5,157.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-0' d='M 145.5,21.3 L 190.9,100.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 145.5,178.7 L 54.5,178.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 54.5,178.7 L 9.1,100.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 63.5,157.8 L 31.7,102.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 9.1,100.0 L 54.5,21.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 54.5,21.3 L 145.5,21.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 68.2,39.5 L 131.8,39.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
</svg>
</div>




| molecule | AMSR |
| ---  | ---  |
| [acetaldehyde](https://en.wikipedia.org/wiki/Acetaldehyde) - only one double bond added | cco |

<div><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='200px' height='200px' viewBox='0 0 200 200'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='200.0' height='200.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 9.1,117.8 L 94.1,68.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 89.2,77.3 L 124.1,97.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 124.1,97.4 L 159.0,117.5' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 99.0,60.3 L 133.9,80.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 133.9,80.4 L 168.8,100.5' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-2' d='M 166.3 117.9
Q 166.3 111.2, 169.6 107.5
Q 172.9 103.8, 179.1 103.8
Q 185.3 103.8, 188.6 107.5
Q 191.9 111.2, 191.9 117.9
Q 191.9 124.7, 188.5 128.5
Q 185.2 132.3, 179.1 132.3
Q 173.0 132.3, 169.6 128.5
Q 166.3 124.7, 166.3 117.9
M 179.1 129.2
Q 183.3 129.2, 185.6 126.3
Q 187.9 123.5, 187.9 117.9
Q 187.9 112.5, 185.6 109.7
Q 183.3 106.9, 179.1 106.9
Q 174.9 106.9, 172.5 109.7
Q 170.3 112.4, 170.3 117.9
Q 170.3 123.5, 172.5 126.3
Q 174.9 129.2, 179.1 129.2
' fill='#FF0000'/>
</svg>
</div>


Note that an oxygen with two neighbors or a nitrogen with three in an aromatic ring
is still denoted by a capital (not a lowercase) symbol,
although sp<sup>2</sup>-hybridized,
since its coordination number is still equal to its valence.


| molecule | AMSR |
| ---  | ---  |
| [furan](https://en.wikipedia.org/wiki/Furan) | ccccO5 |

<div><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='200px' height='200px' viewBox='0 0 200 200'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='200.0' height='200.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 180.1,93.3 L 118.9,9.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 154.1,92.9 L 111.2,34.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-0' d='M 131.2,160.6 L 155.6,126.9' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-0' d='M 155.6,126.9 L 180.1,93.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 118.9,9.1 L 19.9,41.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 19.9,41.3 L 19.9,145.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 40.7,56.9 L 40.7,129.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 19.9,145.3 L 61.6,158.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 61.6,158.9 L 103.3,172.4' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-4' d='M 105.9 177.6
Q 105.9 170.8, 109.3 167.0
Q 112.6 163.2, 118.9 163.2
Q 125.2 163.2, 128.5 167.0
Q 131.9 170.8, 131.9 177.6
Q 131.9 184.5, 128.5 188.4
Q 125.1 192.3, 118.9 192.3
Q 112.7 192.3, 109.3 188.4
Q 105.9 184.5, 105.9 177.6
M 118.9 189.1
Q 123.2 189.1, 125.5 186.2
Q 127.9 183.3, 127.9 177.6
Q 127.9 172.0, 125.5 169.2
Q 123.2 166.4, 118.9 166.4
Q 114.6 166.4, 112.2 169.2
Q 109.9 172.0, 109.9 177.6
Q 109.9 183.3, 112.2 186.2
Q 114.6 189.1, 118.9 189.1
' fill='#FF0000'/>
</svg>
</div>




| molecule | AMSR |
| ---  | ---  |
| [pyrrole](https://en.wikipedia.org/wiki/Pyrrole) | ccccN5 |

<div><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='200px' height='200px' viewBox='0 0 200 200'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='200.0' height='200.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 167.5,80.1 L 115.9,9.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 145.6,79.8 L 109.5,30.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-0' d='M 126.6,136.5 L 147.1,108.3' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-0' d='M 147.1,108.3 L 167.5,80.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 115.9,9.1 L 32.5,36.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 32.5,36.2 L 32.5,124.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 50.0,49.4 L 50.0,110.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 32.5,124.0 L 68.4,135.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 68.4,135.7 L 104.3,147.4' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-4' d='M 110.4 138.7
L 118.6 151.9
Q 119.4 153.2, 120.7 155.5
Q 122.0 157.9, 122.1 158.0
L 122.1 138.7
L 125.4 138.7
L 125.4 163.6
L 122.0 163.6
L 113.2 149.2
Q 112.2 147.5, 111.1 145.5
Q 110.1 143.6, 109.7 143.0
L 109.7 163.6
L 106.5 163.6
L 106.5 138.7
L 110.4 138.7
' fill='#0000FF'/>
<path  class='atom-4' d='M 106.2 166.0
L 109.6 166.0
L 109.6 176.6
L 122.3 176.6
L 122.3 166.0
L 125.7 166.0
L 125.7 190.9
L 122.3 190.9
L 122.3 179.4
L 109.6 179.4
L 109.6 190.9
L 106.2 190.9
L 106.2 166.0
' fill='#0000FF'/>
</svg>
</div>



### Triple bonds (sp centers)
Atoms with two fewer neighbors than their valence are designated by a trailing colon (`:`)
can make triple bonds (or more than one double bond).


| molecule | AMSR |
| ---  | ---  |
| [hydrogen cyanide](https://en.wikipedia.org/wiki/Hydrogen_cyanide) | C\:N\: |

<div><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='200px' height='200px' viewBox='0 0 200 200'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='200.0' height='200.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 173.1,100.0 L 112.0,100.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 112.0,100.0 L 50.8,100.0' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 173.1,72.7 L 112.0,72.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 112.0,72.7 L 50.8,72.7' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 173.1,127.3 L 112.0,127.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 112.0,127.3 L 50.8,127.3' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-1' d='M 30.4 85.8
L 39.7 100.8
Q 40.6 102.3, 42.1 105.0
Q 43.6 107.7, 43.7 107.8
L 43.7 85.8
L 47.4 85.8
L 47.4 114.2
L 43.6 114.2
L 33.6 97.8
Q 32.4 95.8, 31.2 93.6
Q 30.0 91.4, 29.6 90.8
L 29.6 114.2
L 26.0 114.2
L 26.0 85.8
L 30.4 85.8
' fill='#0000FF'/>
</svg>
</div>




| molecule | AMSR |
| ---  | ---  |
| [carbon dioxide](https://en.wikipedia.org/wiki/Carbon_dioxide) | oC\:o |

<div><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='200px' height='200px' viewBox='0 0 200 200'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='200.0' height='200.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 168.0,92.0 L 139.2,92.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 139.2,92.0 L 110.5,92.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 168.0,108.0 L 139.2,108.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 139.2,108.0 L 110.5,108.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 89.2,92.0 L 60.5,92.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 60.5,92.0 L 31.7,92.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 89.2,108.0 L 60.5,108.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 60.5,108.0 L 31.7,108.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-0' d='M 170.0 100.1
Q 170.0 94.6, 172.7 91.5
Q 175.4 88.5, 180.5 88.5
Q 185.5 88.5, 188.2 91.5
Q 190.9 94.6, 190.9 100.1
Q 190.9 105.6, 188.2 108.8
Q 185.4 111.9, 180.5 111.9
Q 175.4 111.9, 172.7 108.8
Q 170.0 105.6, 170.0 100.1
M 180.5 109.3
Q 183.9 109.3, 185.8 107.0
Q 187.7 104.6, 187.7 100.1
Q 187.7 95.6, 185.8 93.3
Q 183.9 91.1, 180.5 91.1
Q 177.0 91.1, 175.1 93.3
Q 173.2 95.6, 173.2 100.1
Q 173.2 104.7, 175.1 107.0
Q 177.0 109.3, 180.5 109.3
' fill='#FF0000'/>
<path  class='atom-1' d='M 91.2 100.1
Q 91.2 94.4, 93.9 91.5
Q 96.5 88.5, 101.6 88.5
Q 106.3 88.5, 108.8 91.8
L 106.7 93.5
Q 104.8 91.1, 101.6 91.1
Q 98.1 91.1, 96.3 93.4
Q 94.5 95.7, 94.5 100.1
Q 94.5 104.6, 96.4 106.9
Q 98.3 109.2, 101.9 109.2
Q 104.4 109.2, 107.4 107.7
L 108.3 110.1
Q 107.1 110.9, 105.3 111.4
Q 103.5 111.8, 101.5 111.8
Q 96.5 111.8, 93.9 108.8
Q 91.2 105.8, 91.2 100.1
' fill='#000000'/>
<path  class='atom-2' d='M 9.1 100.1
Q 9.1 94.6, 11.8 91.5
Q 14.5 88.5, 19.5 88.5
Q 24.6 88.5, 27.3 91.5
Q 30.0 94.6, 30.0 100.1
Q 30.0 105.6, 27.3 108.8
Q 24.5 111.9, 19.5 111.9
Q 14.5 111.9, 11.8 108.8
Q 9.1 105.6, 9.1 100.1
M 19.5 109.3
Q 23.0 109.3, 24.9 107.0
Q 26.8 104.6, 26.8 100.1
Q 26.8 95.6, 24.9 93.3
Q 23.0 91.1, 19.5 91.1
Q 16.1 91.1, 14.2 93.3
Q 12.3 95.6, 12.3 100.1
Q 12.3 104.7, 14.2 107.0
Q 16.1 109.3, 19.5 109.3
' fill='#FF0000'/>
</svg>
</div>



### Atoms with multiple valences
Atoms denoted by their symbol alone are assumed to have their lowest possible valence
(for instance, two for sulfur).  Higher valences are denoted by one or more exclamation points (`!`).


## Developing
This repo uses pre-commit, so after cloning run `pip install -r requirements.txt` and `pre-commit install` prior to committing.
