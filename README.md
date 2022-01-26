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


| AMSR | molecule |
| --- | --- |
C | [methane](https://en.wikipedia.org/wiki/Methane)
O | [water](https://en.wikipedia.org/wiki/Water)
[Cl] | [hydrochloric acid](https://en.wikipedia.org/wiki/Hydrochloric_acid)


<div>
<span style="margin:20px"><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='160px' height='160px' viewBox='0 0 160 160'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='160.0' height='160.0' x='0.0' y='0.0'> </rect>
<path  class='atom-0' d='M 69.1 80.1
Q 69.1 73.1, 72.4 69.4
Q 75.7 65.7, 82.0 65.7
Q 87.8 65.7, 90.9 69.8
L 88.3 72.0
Q 86.0 69.0, 82.0 69.0
Q 77.7 69.0, 75.4 71.8
Q 73.2 74.7, 73.2 80.1
Q 73.2 85.7, 75.5 88.6
Q 77.8 91.5, 82.4 91.5
Q 85.5 91.5, 89.2 89.6
L 90.3 92.6
Q 88.8 93.6, 86.6 94.1
Q 84.3 94.7, 81.8 94.7
Q 75.7 94.7, 72.4 90.9
Q 69.1 87.2, 69.1 80.1
' fill='#000000'/>
<path  class='atom-0' d='M 94.3 66.0
L 98.2 66.0
L 98.2 78.0
L 112.6 78.0
L 112.6 66.0
L 116.5 66.0
L 116.5 94.3
L 112.6 94.3
L 112.6 81.2
L 98.2 81.2
L 98.2 94.3
L 94.3 94.3
L 94.3 66.0
' fill='#000000'/>
<path  class='atom-0' d='M 129.4 102.1
L 131.7 102.1
L 131.7 104.2
L 129.4 104.2
L 129.4 108.5
L 126.9 108.5
L 126.9 104.2
L 117.3 104.2
L 117.3 102.5
L 125.5 89.8
L 129.4 89.8
L 129.4 102.1
M 120.3 102.1
L 126.9 102.1
L 126.9 91.5
L 120.3 102.1
' fill='#000000'/>
</svg>
</span>
<span style="margin:20px"><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='160px' height='160px' viewBox='0 0 160 160'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='160.0' height='160.0' x='0.0' y='0.0'> </rect>
<path  class='atom-0' d='M 28.9 66.0
L 32.8 66.0
L 32.8 78.0
L 47.2 78.0
L 47.2 66.0
L 51.1 66.0
L 51.1 94.3
L 47.2 94.3
L 47.2 81.2
L 32.8 81.2
L 32.8 94.3
L 28.9 94.3
L 28.9 66.0
' fill='#FF0000'/>
<path  class='atom-0' d='M 52.5 93.3
Q 53.1 91.6, 54.8 90.6
Q 56.4 89.6, 58.7 89.6
Q 61.5 89.6, 63.1 91.1
Q 64.7 92.6, 64.7 95.4
Q 64.7 98.1, 62.6 100.7
Q 60.6 103.3, 56.4 106.4
L 65.0 106.4
L 65.0 108.5
L 52.4 108.5
L 52.4 106.7
Q 55.9 104.2, 58.0 102.4
Q 60.0 100.5, 61.0 98.9
Q 62.0 97.2, 62.0 95.5
Q 62.0 93.7, 61.1 92.7
Q 60.2 91.7, 58.7 91.7
Q 57.2 91.7, 56.2 92.3
Q 55.2 92.9, 54.5 94.3
L 52.5 93.3
' fill='#FF0000'/>
<path  class='atom-0' d='M 67.0 80.1
Q 67.0 73.3, 70.4 69.5
Q 73.7 65.7, 80.0 65.7
Q 86.3 65.7, 89.6 69.5
Q 93.0 73.3, 93.0 80.1
Q 93.0 87.0, 89.6 90.9
Q 86.2 94.8, 80.0 94.8
Q 73.8 94.8, 70.4 90.9
Q 67.0 87.0, 67.0 80.1
M 80.0 91.6
Q 84.3 91.6, 86.6 88.7
Q 89.0 85.8, 89.0 80.1
Q 89.0 74.5, 86.6 71.7
Q 84.3 68.9, 80.0 68.9
Q 75.7 68.9, 73.3 71.7
Q 71.0 74.5, 71.0 80.1
Q 71.0 85.8, 73.3 88.7
Q 75.7 91.6, 80.0 91.6
' fill='#FF0000'/>
</svg>
</span>
<span style="margin:20px"><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='160px' height='160px' viewBox='0 0 160 160'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='160.0' height='160.0' x='0.0' y='0.0'> </rect>
<path  class='atom-0' d='M 28.7 66.9
L 32.6 66.9
L 32.6 78.9
L 47.0 78.9
L 47.0 66.9
L 50.9 66.9
L 50.9 95.2
L 47.0 95.2
L 47.0 82.1
L 32.6 82.1
L 32.6 95.2
L 28.7 95.2
L 28.7 66.9
' fill='#00CC00'/>
<path  class='atom-0' d='M 52.9 81.0
Q 52.9 73.9, 56.2 70.3
Q 59.5 66.5, 65.8 66.5
Q 71.6 66.5, 74.7 70.7
L 72.1 72.8
Q 69.8 69.8, 65.8 69.8
Q 61.5 69.8, 59.2 72.7
Q 57.0 75.5, 57.0 81.0
Q 57.0 86.6, 59.3 89.5
Q 61.6 92.3, 66.2 92.3
Q 69.3 92.3, 73.0 90.5
L 74.1 93.5
Q 72.6 94.4, 70.4 95.0
Q 68.1 95.5, 65.6 95.5
Q 59.5 95.5, 56.2 91.8
Q 52.9 88.0, 52.9 81.0
' fill='#00CC00'/>
<path  class='atom-0' d='M 78.2 64.8
L 81.8 64.8
L 81.8 95.2
L 78.2 95.2
L 78.2 64.8
' fill='#00CC00'/>
</svg>
</span>
</div>


### Chains
Each atom in a chain is bonded to the most recently added atom that
can still make bonds, according to its valence. Hydrogens may be added
explicitly like any other atom.  In the following
example, the fluorines are added to the second carbon.  The chlorine
is then added to the first carbon, since the second can no longer bond.

| AMSR | molecule |
| --- | --- |
CCFFF[Cl] | [2-chloro-1,1,1-trifluoroethane](https://pubchem.ncbi.nlm.nih.gov/compound/2-Chloro-1_1_1-trifluoroethane)


<div>
<span style="margin:20px"><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='160px' height='160px' viewBox='0 0 160 160'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='160.0' height='160.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 63.3,104.6 L 105.9,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-0 atom-5' d='M 63.3,104.6 L 43.0,92.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-0 atom-5' d='M 43.0,92.9 L 22.8,81.2' style='fill:none;fill-rule:evenodd;stroke:#00CC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 105.9,80.0 L 124.6,69.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 124.6,69.2 L 143.2,58.5' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-1 atom-3' d='M 105.9,80.0 L 115.9,97.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-1 atom-3' d='M 115.9,97.2 L 125.8,114.4' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-1 atom-4' d='M 105.9,80.0 L 96.0,62.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-1 atom-4' d='M 96.0,62.8 L 86.0,45.6' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-2' d='M 144.4 48.4
L 152.7 48.4
L 152.7 50.0
L 146.3 50.0
L 146.3 54.2
L 152.0 54.2
L 152.0 55.8
L 146.3 55.8
L 146.3 62.3
L 144.4 62.3
L 144.4 48.4
' fill='#33CCCC'/>
<path  class='atom-3' d='M 126.4 115.7
L 134.7 115.7
L 134.7 117.3
L 128.3 117.3
L 128.3 121.5
L 134.0 121.5
L 134.0 123.1
L 128.3 123.1
L 128.3 129.6
L 126.4 129.6
L 126.4 115.7
' fill='#33CCCC'/>
<path  class='atom-4' d='M 77.2 30.4
L 85.5 30.4
L 85.5 32.0
L 79.0 32.0
L 79.0 36.2
L 84.7 36.2
L 84.7 37.8
L 79.0 37.8
L 79.0 44.3
L 77.2 44.3
L 77.2 30.4
' fill='#33CCCC'/>
<path  class='atom-5' d='M 7.3 80.5
Q 7.3 77.0, 8.9 75.2
Q 10.5 73.4, 13.6 73.4
Q 16.5 73.4, 18.0 75.4
L 16.7 76.5
Q 15.6 75.0, 13.6 75.0
Q 11.5 75.0, 10.4 76.4
Q 9.3 77.8, 9.3 80.5
Q 9.3 83.2, 10.4 84.7
Q 11.6 86.1, 13.8 86.1
Q 15.4 86.1, 17.2 85.2
L 17.7 86.6
Q 17.0 87.1, 15.9 87.4
Q 14.8 87.7, 13.6 87.7
Q 10.5 87.7, 8.9 85.8
Q 7.3 83.9, 7.3 80.5
' fill='#00CC00'/>
<path  class='atom-5' d='M 19.7 72.5
L 21.5 72.5
L 21.5 87.5
L 19.7 87.5
L 19.7 72.5
' fill='#00CC00'/>
</svg>
</span>
</div>

If no atoms are available to bond, a new molecule is formed.

| AMSR | molecule |
| --- | --- |
CFFFFCO | [carbon tetrafluoride](https://en.wikipedia.org/wiki/Carbon_tetrafluoride) and [methanol](https://en.wikipedia.org/wiki/Methanol)


<div>
<span style="margin:20px"><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='160px' height='160px' viewBox='0 0 160 160'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='160.0' height='160.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 80.0,63.4 L 101.9,63.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 101.9,63.4 L 123.8,63.4' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 80.0,63.4 L 58.1,63.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 58.1,63.4 L 36.2,63.4' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 80.0,63.4 L 80.0,42.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 80.0,42.9 L 80.0,22.4' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-0 atom-4' d='M 80.0,63.4 L 80.0,84.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-0 atom-4' d='M 80.0,84.0 L 80.0,104.5' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-5 atom-6' d='M 55.4,145.5 L 76.2,145.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-5 atom-6' d='M 76.2,145.5 L 97.0,145.5' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-1' d='M 125.1 56.5
L 133.4 56.5
L 133.4 58.1
L 126.9 58.1
L 126.9 62.3
L 132.6 62.3
L 132.6 63.9
L 126.9 63.9
L 126.9 70.4
L 125.1 70.4
L 125.1 56.5
' fill='#33CCCC'/>
<path  class='atom-2' d='M 26.6 56.5
L 34.9 56.5
L 34.9 58.1
L 28.5 58.1
L 28.5 62.3
L 34.2 62.3
L 34.2 63.9
L 28.5 63.9
L 28.5 70.4
L 26.6 70.4
L 26.6 56.5
' fill='#33CCCC'/>
<path  class='atom-3' d='M 75.9 7.3
L 84.1 7.3
L 84.1 8.9
L 77.7 8.9
L 77.7 13.1
L 83.4 13.1
L 83.4 14.7
L 77.7 14.7
L 77.7 21.2
L 75.9 21.2
L 75.9 7.3
' fill='#33CCCC'/>
<path  class='atom-4' d='M 75.9 105.7
L 84.1 105.7
L 84.1 107.3
L 77.7 107.3
L 77.7 111.5
L 83.4 111.5
L 83.4 113.1
L 77.7 113.1
L 77.7 119.6
L 75.9 119.6
L 75.9 105.7
' fill='#33CCCC'/>
<path  class='atom-6' d='M 98.2 145.5
Q 98.2 142.2, 99.9 140.3
Q 101.5 138.4, 104.6 138.4
Q 107.7 138.4, 109.3 140.3
Q 111.0 142.2, 111.0 145.5
Q 111.0 148.9, 109.3 150.8
Q 107.7 152.7, 104.6 152.7
Q 101.5 152.7, 99.9 150.8
Q 98.2 148.9, 98.2 145.5
M 104.6 151.2
Q 106.7 151.2, 107.9 149.7
Q 109.0 148.3, 109.0 145.5
Q 109.0 142.8, 107.9 141.4
Q 106.7 140.0, 104.6 140.0
Q 102.5 140.0, 101.3 141.4
Q 100.2 142.7, 100.2 145.5
Q 100.2 148.3, 101.3 149.7
Q 102.5 151.2, 104.6 151.2
' fill='#FF0000'/>
<path  class='atom-6' d='M 112.7 138.6
L 114.6 138.6
L 114.6 144.5
L 121.7 144.5
L 121.7 138.6
L 123.6 138.6
L 123.6 152.5
L 121.7 152.5
L 121.7 146.1
L 114.6 146.1
L 114.6 152.5
L 112.7 152.5
L 112.7 138.6
' fill='#FF0000'/>
</svg>
</span>
</div>


### Branches
Branches are formed automatically when atoms can no longer
make bonds.  They can also be made by "capping" or
"saturating" an atom with hydrogens, using a period `.`.
The capping hydrogens are applied to the most
recently-added atom which can still make bonds.
New atoms will then be added to those added earlier, forming a branch.

| AMSR | molecule |
| --- | --- |
CCC.C | [isobutane](https://en.wikipedia.org/wiki/Isobutane)
CC.CC.C.C | [2,2-dimethylbutane](https://pubchem.ncbi.nlm.nih.gov/compound/6403)


<div>
<span style="margin:20px"><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='160px' height='160px' viewBox='0 0 160 160'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='160.0' height='160.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 7.3,143.0 L 80.0,101.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 80.0,101.0 L 152.7,143.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-1 atom-3' d='M 80.0,101.0 L 80.0,17.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
</svg>
</span>
<span style="margin:20px"><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='160px' height='160px' viewBox='0 0 160 160'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='160.0' height='160.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 55.8,52.0 L 7.3,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 55.8,52.0 L 104.2,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 104.2,80.0 L 132.2,31.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-2 atom-4' d='M 104.2,80.0 L 76.2,128.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-2 atom-5' d='M 104.2,80.0 L 152.7,108.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
</svg>
</span>
</div>


### Small rings
Small rings are denoted by one of the digits `3,4,5,` or `6` indicating
the size of the ring.
A new bond is formed between atoms *i* and *j* where
*i* and *j* are the most recently-added atoms which
can still make bonds and when bonded will form a ring of that size.

| AMSR | molecule |
| --- | --- |
CCO3 | [oxirane](https://en.wikipedia.org/wiki/Ethylene_oxide)
CCCCCC6 | [cyclohexane](https://en.wikipedia.org/wiki/Cyclohexane)


<div>
<span style="margin:20px"><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='160px' height='160px' viewBox='0 0 160 160'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='160.0' height='160.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 143.1,73.3 L 28.7,7.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-0' d='M 45.0,129.9 L 94.1,101.6' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-0' d='M 94.1,101.6 L 143.1,73.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 28.7,7.3 L 28.7,64.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 28.7,64.5 L 28.7,121.7' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-2' d='M 15.7 139.4
Q 15.7 132.6, 19.1 128.8
Q 22.4 125.0, 28.7 125.0
Q 35.0 125.0, 38.4 128.8
Q 41.7 132.6, 41.7 139.4
Q 41.7 146.3, 38.3 150.2
Q 34.9 154.1, 28.7 154.1
Q 22.5 154.1, 19.1 150.2
Q 15.7 146.3, 15.7 139.4
M 28.7 150.9
Q 33.0 150.9, 35.4 148.0
Q 37.7 145.1, 37.7 139.4
Q 37.7 133.8, 35.4 131.0
Q 33.0 128.2, 28.7 128.2
Q 24.4 128.2, 22.0 131.0
Q 19.7 133.8, 19.7 139.4
Q 19.7 145.1, 22.0 148.0
Q 24.4 150.9, 28.7 150.9
' fill='#FF0000'/>
</svg>
</span>
<span style="margin:20px"><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='160px' height='160px' viewBox='0 0 160 160'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='160.0' height='160.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 152.7,80.0 L 116.4,143.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-0' d='M 116.4,17.0 L 152.7,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 116.4,143.0 L 43.6,143.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 43.6,143.0 L 7.3,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 7.3,80.0 L 43.6,17.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 43.6,17.0 L 116.4,17.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
</svg>
</span>
</div>


### Large rings
Rings with seven or more members are denoted by a
sequence of digits `3,4,5,` or `6` which are added to give the size of the ring.
If any characters appear between digits, multiple rings are formed
rather than a single large ring.  Whitespace may be used if this is desired
(without adding any additional atoms).

| AMSR | molecule |
| --- | --- |
CCCCCCC43 | [cycloheptane](https://en.wikipedia.org/wiki/Cycloheptane)
CCCCCCCCCCCC66 | [cyclododecane](https://en.wikipedia.org/wiki/Cyclododecane)
CCCCCCCCCCCC6 6 | [1-ethyldecalin](https://pubchem.ncbi.nlm.nih.gov/compound/33053)


<div>
<span style="margin:20px"><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='160px' height='160px' viewBox='0 0 160 160'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='160.0' height='160.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 150.9,80.0 L 122.8,21.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-6 atom-0' d='M 122.8,138.3 L 150.9,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 122.8,21.7 L 59.7,7.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 59.7,7.3 L 9.1,47.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 9.1,47.6 L 9.1,112.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 9.1,112.4 L 59.7,152.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-6' d='M 59.7,152.7 L 122.8,138.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
</svg>
</span>
<span style="margin:20px"><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='160px' height='160px' viewBox='0 0 160 160'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='160.0' height='160.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 152.7,80.0 L 143.0,116.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-11 atom-0' d='M 143.0,43.6 L 152.7,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 143.0,116.4 L 116.4,143.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 116.4,143.0 L 80.0,152.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 80.0,152.7 L 43.6,143.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 43.6,143.0 L 17.0,116.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-6' d='M 17.0,116.4 L 7.3,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-6 atom-7' d='M 7.3,80.0 L 17.0,43.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-7 atom-8' d='M 17.0,43.6 L 43.6,17.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-8 atom-9' d='M 43.6,17.0 L 80.0,7.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-9 atom-10' d='M 80.0,7.3 L 116.4,17.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-10 atom-11' d='M 116.4,17.0 L 143.0,43.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
</svg>
</span>
<span style="margin:20px"><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='160px' height='160px' viewBox='0 0 160 160'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='160.0' height='160.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 7.3,44.6 L 39.2,42.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 39.2,42.4 L 57.0,69.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 57.0,69.0 L 43.0,97.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-12 atom-11 atom-2' d='M 88.9,66.7 L 57.0,69.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 43.0,97.7 L 60.9,124.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 60.9,124.2 L 92.8,122.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-6' d='M 92.8,122.0 L 106.8,93.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-6 atom-7' d='M 106.8,93.3 L 138.7,91.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-11 atom-6' d='M 88.9,66.7 L 106.8,93.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-7 atom-8' d='M 138.7,91.0 L 152.7,62.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-8 atom-9' d='M 152.7,62.3 L 134.9,35.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-9 atom-10' d='M 134.9,35.8 L 103.0,38.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-10 atom-11' d='M 103.0,38.0 L 88.9,66.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
</svg>
</span>
</div>


### Double bonds (sp<sup>2</sup> centers)
Atoms making a single double bond are indicated by changing
the symbol to lowercase (note that lowercase
does not mean "aromatic"; merely, "atom having one
fewer neighbor than its valence.")  Double bonds are
assigned by a matching algorithm.  If a perfect
matching cannot be found (for instance, in the case
of an odd number of contiguous lowercase symbols) a
maximal matching is chosen, non-matched atoms
remain singly bonded, and hydrogens are added.

| AMSR | molecule |
| --- | --- |
co | [formaldehyde](https://en.wikipedia.org/wiki/Formaldehyde)
cccccc6 | [benzene](https://en.wikipedia.org/wiki/Benzene)
cco | [acetaldehyde](https://en.wikipedia.org/wiki/Acetaldehyde) (only one double bond added)


<div>
<span style="margin:20px"><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='160px' height='160px' viewBox='0 0 160 160'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='160.0' height='160.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 19.5,90.9 L 66.2,90.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 66.2,90.9 L 112.9,90.9' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 19.5,69.1 L 66.2,69.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 66.2,69.1 L 112.9,69.1' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-1' d='M 115.6 80.1
Q 115.6 73.3, 119.0 69.5
Q 122.4 65.7, 128.6 65.7
Q 134.9 65.7, 138.3 69.5
Q 141.6 73.3, 141.6 80.1
Q 141.6 87.0, 138.2 90.9
Q 134.8 94.8, 128.6 94.8
Q 122.4 94.8, 119.0 90.9
Q 115.6 87.0, 115.6 80.1
M 128.6 91.6
Q 133.0 91.6, 135.3 88.7
Q 137.6 85.8, 137.6 80.1
Q 137.6 74.5, 135.3 71.7
Q 133.0 68.9, 128.6 68.9
Q 124.3 68.9, 122.0 71.7
Q 119.6 74.5, 119.6 80.1
Q 119.6 85.8, 122.0 88.7
Q 124.3 91.6, 128.6 91.6
' fill='#FF0000'/>
</svg>
</span>
<span style="margin:20px"><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='160px' height='160px' viewBox='0 0 160 160'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='160.0' height='160.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 152.7,80.0 L 116.4,143.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 134.7,82.2 L 109.2,126.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-0' d='M 116.4,17.0 L 152.7,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 116.4,143.0 L 43.6,143.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 43.6,143.0 L 7.3,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 50.8,126.3 L 25.3,82.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 7.3,80.0 L 43.6,17.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 43.6,17.0 L 116.4,17.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 54.5,31.6 L 105.5,31.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
</svg>
</span>
<span style="margin:20px"><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='160px' height='160px' viewBox='0 0 160 160'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='160.0' height='160.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 7.3,93.8 L 74.9,54.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 71.0,61.5 L 98.8,77.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 98.8,77.5 L 126.6,93.5' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 78.8,47.9 L 106.6,64.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 106.6,64.0 L 134.4,80.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-2' d='M 132.4 93.8
Q 132.4 88.5, 135.0 85.5
Q 137.7 82.6, 142.6 82.6
Q 147.5 82.6, 150.1 85.5
Q 152.7 88.5, 152.7 93.8
Q 152.7 99.2, 150.1 102.3
Q 147.4 105.3, 142.6 105.3
Q 137.7 105.3, 135.0 102.3
Q 132.4 99.2, 132.4 93.8
M 142.6 102.8
Q 145.9 102.8, 147.8 100.5
Q 149.6 98.3, 149.6 93.8
Q 149.6 89.5, 147.8 87.3
Q 145.9 85.1, 142.6 85.1
Q 139.2 85.1, 137.4 87.3
Q 135.5 89.5, 135.5 93.8
Q 135.5 98.3, 137.4 100.5
Q 139.2 102.8, 142.6 102.8
' fill='#FF0000'/>
</svg>
</span>
</div>


Note that an oxygen with two neighbors or a nitrogen with three in an aromatic ring
is still denoted by a capital (not a lowercase) symbol,
although sp<sup>2</sup>-hybridized,
since its coordination number is still equal to its valence.

| AMSR | molecule |
| --- | --- |
ccccO5 | [furan](https://en.wikipedia.org/wiki/Furan)
ccccN5 | [pyrrole](https://en.wikipedia.org/wiki/Pyrrole)


<div>
<span style="margin:20px"><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='160px' height='160px' viewBox='0 0 160 160'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='160.0' height='160.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 143.4,73.9 L 95.0,7.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 122.8,73.6 L 88.9,27.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-0' d='M 105.0,126.7 L 124.2,100.3' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-0' d='M 124.2,100.3 L 143.4,73.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 95.0,7.3 L 16.6,32.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 16.6,32.7 L 16.6,115.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 33.1,45.1 L 33.1,102.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 16.6,115.1 L 49.4,125.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 49.4,125.8 L 82.2,136.4' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-4' d='M 84.3 140.6
Q 84.3 135.0, 87.0 131.9
Q 89.8 128.8, 95.0 128.8
Q 100.1 128.8, 102.9 131.9
Q 105.7 135.0, 105.7 140.6
Q 105.7 146.3, 102.9 149.5
Q 100.1 152.7, 95.0 152.7
Q 89.8 152.7, 87.0 149.5
Q 84.3 146.3, 84.3 140.6
M 95.0 150.1
Q 98.5 150.1, 100.4 147.7
Q 102.4 145.3, 102.4 140.6
Q 102.4 136.1, 100.4 133.7
Q 98.5 131.4, 95.0 131.4
Q 91.4 131.4, 89.5 133.7
Q 87.5 136.0, 87.5 140.6
Q 87.5 145.3, 89.5 147.7
Q 91.4 150.1, 95.0 150.1
' fill='#FF0000'/>
</svg>
</span>
<span style="margin:20px"><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='160px' height='160px' viewBox='0 0 160 160'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='160.0' height='160.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 134.0,64.1 L 92.8,7.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 116.5,63.8 L 87.6,24.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-0' d='M 101.3,109.2 L 117.6,86.6' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-0' d='M 117.6,86.6 L 134.0,64.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 92.8,7.3 L 26.0,29.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 26.0,29.0 L 26.0,99.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 40.0,39.5 L 40.0,88.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 26.0,99.2 L 54.7,108.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 54.7,108.5 L 83.5,117.9' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-4' d='M 88.4 111.0
L 94.9 121.5
Q 95.5 122.5, 96.6 124.4
Q 97.6 126.3, 97.7 126.4
L 97.7 111.0
L 100.3 111.0
L 100.3 130.8
L 97.6 130.8
L 90.6 119.3
Q 89.8 118.0, 88.9 116.4
Q 88.1 114.9, 87.8 114.4
L 87.8 130.8
L 85.2 130.8
L 85.2 111.0
L 88.4 111.0
' fill='#0000FF'/>
<path  class='atom-4' d='M 85.0 132.8
L 87.7 132.8
L 87.7 141.3
L 97.8 141.3
L 97.8 132.8
L 100.5 132.8
L 100.5 152.7
L 97.8 152.7
L 97.8 143.5
L 87.7 143.5
L 87.7 152.7
L 85.0 152.7
L 85.0 132.8
' fill='#0000FF'/>
</svg>
</span>
</div>


### Triple bonds (sp centers)
Atoms with two fewer neighbors than their valence are designated by a trailing colon `:`
can make triple bonds (or more than one double bond).

| AMSR | molecule |
| --- | --- |
C\:N\: | [hydrogen cyanide](https://en.wikipedia.org/wiki/Hydrogen_cyanide)
oC\:o | [carbon dioxide](https://en.wikipedia.org/wiki/Carbon_dioxide)


<div>
<span style="margin:20px"><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='160px' height='160px' viewBox='0 0 160 160'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='160.0' height='160.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 139.4,80.0 L 91.6,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 91.6,80.0 L 43.8,80.0' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 139.4,58.2 L 91.6,58.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 91.6,58.2 L 43.8,58.2' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 139.4,101.8 L 91.6,101.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 91.6,101.8 L 43.8,101.8' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-1' d='M 24.1 65.8
L 33.4 80.8
Q 34.3 82.3, 35.8 85.0
Q 37.2 87.7, 37.3 87.8
L 37.3 65.8
L 41.1 65.8
L 41.1 94.2
L 37.2 94.2
L 27.2 77.8
Q 26.1 75.8, 24.8 73.6
Q 23.6 71.4, 23.3 70.8
L 23.3 94.2
L 19.6 94.2
L 19.6 65.8
L 24.1 65.8
' fill='#0000FF'/>
</svg>
</span>
<span style="margin:20px"><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='160px' height='160px' viewBox='0 0 160 160'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='160.0' height='160.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 134.4,73.6 L 111.4,73.6' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 111.4,73.6 L 88.4,73.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 134.4,86.4 L 111.4,86.4' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 111.4,86.4 L 88.4,86.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 71.4,73.6 L 48.4,73.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 48.4,73.6 L 25.4,73.6' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 71.4,86.4 L 48.4,86.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 48.4,86.4 L 25.4,86.4' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-0' d='M 136.0 80.1
Q 136.0 75.7, 138.2 73.2
Q 140.3 70.8, 144.4 70.8
Q 148.4 70.8, 150.6 73.2
Q 152.7 75.7, 152.7 80.1
Q 152.7 84.5, 150.5 87.0
Q 148.4 89.5, 144.4 89.5
Q 140.3 89.5, 138.2 87.0
Q 136.0 84.5, 136.0 80.1
M 144.4 87.4
Q 147.1 87.4, 148.6 85.6
Q 150.2 83.7, 150.2 80.1
Q 150.2 76.5, 148.6 74.7
Q 147.1 72.8, 144.4 72.8
Q 141.6 72.8, 140.1 74.6
Q 138.6 76.4, 138.6 80.1
Q 138.6 83.7, 140.1 85.6
Q 141.6 87.4, 144.4 87.4
' fill='#FF0000'/>
<path  class='atom-1' d='M 73.0 80.1
Q 73.0 75.5, 75.1 73.2
Q 77.2 70.8, 81.3 70.8
Q 85.0 70.8, 87.0 73.4
L 85.3 74.8
Q 83.9 72.9, 81.3 72.9
Q 78.5 72.9, 77.0 74.7
Q 75.6 76.6, 75.6 80.1
Q 75.6 83.7, 77.1 85.5
Q 78.6 87.4, 81.5 87.4
Q 83.6 87.4, 85.9 86.2
L 86.6 88.1
Q 85.7 88.7, 84.2 89.1
Q 82.8 89.4, 81.2 89.4
Q 77.2 89.4, 75.1 87.0
Q 73.0 84.6, 73.0 80.1
' fill='#000000'/>
<path  class='atom-2' d='M 7.3 80.1
Q 7.3 75.7, 9.4 73.2
Q 11.6 70.8, 15.6 70.8
Q 19.7 70.8, 21.8 73.2
Q 24.0 75.7, 24.0 80.1
Q 24.0 84.5, 21.8 87.0
Q 19.6 89.5, 15.6 89.5
Q 11.6 89.5, 9.4 87.0
Q 7.3 84.5, 7.3 80.1
M 15.6 87.4
Q 18.4 87.4, 19.9 85.6
Q 21.4 83.7, 21.4 80.1
Q 21.4 76.5, 19.9 74.7
Q 18.4 72.8, 15.6 72.8
Q 12.9 72.8, 11.3 74.6
Q 9.8 76.4, 9.8 80.1
Q 9.8 83.7, 11.3 85.6
Q 12.9 87.4, 15.6 87.4
' fill='#FF0000'/>
</svg>
</span>
</div>


### Hypervalent atoms
Atoms denoted by their symbol alone are assumed to have their lowest possible valence
(for instance, two for sulfur).  Higher valences are denoted by one or more exclamation points `!`.

| AMSR | molecule |
| --- | --- |
CSC | [dimethyl sulfide](https://en.wikipedia.org/wiki/Dimethyl_sulfide)
Cs!oC | [dimethyl sulfoxide](https://en.wikipedia.org/wiki/Dimethyl_sulfoxide)
S!!FFFFFF | [sulfur hexafluoride](https://en.wikipedia.org/wiki/Sulfur_hexafluoride)


<div>
<span style="margin:20px"><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='160px' height='160px' viewBox='0 0 160 160'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='160.0' height='160.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 7.3,106.9 L 38.8,88.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 38.8,88.7 L 70.3,70.5' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 89.7,70.5 L 121.2,88.7' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 121.2,88.7 L 152.7,106.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-1' d='M 73.3 73.1
Q 73.6 73.2, 74.7 73.7
Q 75.8 74.1, 77.0 74.4
Q 78.2 74.7, 79.4 74.7
Q 81.7 74.7, 83.0 73.6
Q 84.3 72.5, 84.3 70.6
Q 84.3 69.3, 83.6 68.5
Q 83.0 67.7, 82.0 67.3
Q 81.0 66.8, 79.3 66.3
Q 77.2 65.7, 75.9 65.1
Q 74.7 64.5, 73.8 63.2
Q 72.9 61.9, 72.9 59.8
Q 72.9 56.8, 74.9 54.9
Q 76.9 53.1, 81.0 53.1
Q 83.7 53.1, 86.9 54.4
L 86.1 57.0
Q 83.2 55.8, 81.1 55.8
Q 78.8 55.8, 77.5 56.8
Q 76.2 57.7, 76.2 59.4
Q 76.2 60.6, 76.9 61.4
Q 77.5 62.2, 78.5 62.6
Q 79.5 63.1, 81.1 63.6
Q 83.2 64.2, 84.5 64.9
Q 85.8 65.6, 86.7 66.9
Q 87.6 68.3, 87.6 70.6
Q 87.6 73.9, 85.4 75.7
Q 83.2 77.4, 79.6 77.4
Q 77.4 77.4, 75.8 77.0
Q 74.3 76.5, 72.4 75.8
L 73.3 73.1
' fill='#CCCC00'/>
</svg>
</span>
<span style="margin:20px"><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='160px' height='160px' viewBox='0 0 160 160'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='160.0' height='160.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 7.3,149.0 L 38.8,130.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 38.8,130.8 L 70.3,112.6' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 88.4,93.1 L 88.4,65.1' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 88.4,65.1 L 88.4,37.2' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 71.6,93.1 L 71.6,65.1' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 71.6,65.1 L 71.6,37.2' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-1 atom-3' d='M 89.7,112.6 L 121.2,130.8' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-1 atom-3' d='M 121.2,130.8 L 152.7,149.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-1' d='M 73.3 115.2
Q 73.6 115.3, 74.7 115.7
Q 75.8 116.2, 77.0 116.5
Q 78.2 116.8, 79.4 116.8
Q 81.7 116.8, 83.0 115.7
Q 84.3 114.6, 84.3 112.7
Q 84.3 111.4, 83.6 110.6
Q 83.0 109.8, 82.0 109.3
Q 81.0 108.9, 79.3 108.4
Q 77.2 107.7, 75.9 107.1
Q 74.7 106.5, 73.8 105.3
Q 72.9 104.0, 72.9 101.8
Q 72.9 98.8, 74.9 97.0
Q 76.9 95.1, 81.0 95.1
Q 83.7 95.1, 86.9 96.5
L 86.1 99.0
Q 83.2 97.9, 81.1 97.9
Q 78.8 97.9, 77.5 98.8
Q 76.2 99.8, 76.2 101.4
Q 76.2 102.7, 76.9 103.5
Q 77.5 104.3, 78.5 104.7
Q 79.5 105.1, 81.1 105.6
Q 83.2 106.3, 84.5 107.0
Q 85.8 107.6, 86.7 109.0
Q 87.6 110.4, 87.6 112.7
Q 87.6 116.0, 85.4 117.8
Q 83.2 119.5, 79.6 119.5
Q 77.4 119.5, 75.8 119.0
Q 74.3 118.6, 72.4 117.8
L 73.3 115.2
' fill='#CCCC00'/>
<path  class='atom-2' d='M 69.1 23.1
Q 69.1 17.4, 71.9 14.2
Q 74.7 11.0, 80.0 11.0
Q 85.3 11.0, 88.1 14.2
Q 90.9 17.4, 90.9 23.1
Q 90.9 28.9, 88.1 32.2
Q 85.2 35.4, 80.0 35.4
Q 74.8 35.4, 71.9 32.2
Q 69.1 28.9, 69.1 23.1
M 80.0 32.7
Q 83.6 32.7, 85.6 30.3
Q 87.6 27.9, 87.6 23.1
Q 87.6 18.4, 85.6 16.1
Q 83.6 13.7, 80.0 13.7
Q 76.4 13.7, 74.4 16.0
Q 72.4 18.4, 72.4 23.1
Q 72.4 27.9, 74.4 30.3
Q 76.4 32.7, 80.0 32.7
' fill='#FF0000'/>
</svg>
</span>
<span style="margin:20px"><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='160px' height='160px' viewBox='0 0 160 160'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='160.0' height='160.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 72.6,84.3 L 52.1,96.1' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 52.1,96.1 L 31.6,108.0' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 80.0,91.1 L 80.0,112.2' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 80.0,112.2 L 80.0,133.4' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 87.4,84.3 L 107.9,96.1' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 107.9,96.1 L 128.4,108.0' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-0 atom-4' d='M 80.0,69.4 L 80.0,48.0' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-0 atom-4' d='M 80.0,48.0 L 80.0,26.6' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-0 atom-5' d='M 87.4,75.7 L 107.9,63.9' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-0 atom-5' d='M 107.9,63.9 L 128.4,52.0' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-0 atom-6' d='M 72.6,75.7 L 52.1,63.9' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-0 atom-6' d='M 52.1,63.9 L 31.6,52.0' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-0' d='M 74.9 86.2
Q 75.1 86.3, 75.9 86.6
Q 76.8 87.0, 77.7 87.2
Q 78.6 87.4, 79.6 87.4
Q 81.3 87.4, 82.3 86.6
Q 83.3 85.8, 83.3 84.3
Q 83.3 83.3, 82.8 82.7
Q 82.3 82.1, 81.5 81.8
Q 80.7 81.4, 79.5 81.0
Q 77.9 80.6, 76.9 80.1
Q 75.9 79.6, 75.3 78.7
Q 74.6 77.7, 74.6 76.1
Q 74.6 73.8, 76.1 72.4
Q 77.7 71.0, 80.7 71.0
Q 82.8 71.0, 85.2 72.0
L 84.6 74.0
Q 82.4 73.1, 80.8 73.1
Q 79.1 73.1, 78.1 73.8
Q 77.1 74.5, 77.1 75.8
Q 77.1 76.7, 77.6 77.3
Q 78.1 77.9, 78.9 78.2
Q 79.6 78.6, 80.8 79.0
Q 82.4 79.5, 83.4 80.0
Q 84.4 80.5, 85.1 81.5
Q 85.8 82.5, 85.8 84.3
Q 85.8 86.8, 84.1 88.2
Q 82.4 89.5, 79.7 89.5
Q 78.1 89.5, 76.8 89.1
Q 75.6 88.8, 74.2 88.2
L 74.9 86.2
' fill='#CCCC00'/>
<path  class='atom-1' d='M 19.5 102.8
L 30.2 102.8
L 30.2 104.9
L 21.9 104.9
L 21.9 110.4
L 29.3 110.4
L 29.3 112.5
L 21.9 112.5
L 21.9 120.9
L 19.5 120.9
L 19.5 102.8
' fill='#33CCCC'/>
<path  class='atom-2' d='M 74.6 134.7
L 85.4 134.7
L 85.4 136.7
L 77.1 136.7
L 77.1 142.2
L 84.4 142.2
L 84.4 144.3
L 77.1 144.3
L 77.1 152.7
L 74.6 152.7
L 74.6 134.7
' fill='#33CCCC'/>
<path  class='atom-3' d='M 129.8 102.8
L 140.5 102.8
L 140.5 104.9
L 132.2 104.9
L 132.2 110.4
L 139.6 110.4
L 139.6 112.5
L 132.2 112.5
L 132.2 120.9
L 129.8 120.9
L 129.8 102.8
' fill='#33CCCC'/>
<path  class='atom-4' d='M 74.6 7.3
L 85.4 7.3
L 85.4 9.3
L 77.1 9.3
L 77.1 14.8
L 84.4 14.8
L 84.4 16.9
L 77.1 16.9
L 77.1 25.3
L 74.6 25.3
L 74.6 7.3
' fill='#33CCCC'/>
<path  class='atom-5' d='M 129.8 39.1
L 140.5 39.1
L 140.5 41.2
L 132.2 41.2
L 132.2 46.7
L 139.6 46.7
L 139.6 48.8
L 132.2 48.8
L 132.2 57.2
L 129.8 57.2
L 129.8 39.1
' fill='#33CCCC'/>
<path  class='atom-6' d='M 19.5 39.1
L 30.2 39.1
L 30.2 41.2
L 21.9 41.2
L 21.9 46.7
L 29.3 46.7
L 29.3 48.8
L 21.9 48.8
L 21.9 57.2
L 19.5 57.2
L 19.5 39.1
' fill='#33CCCC'/>
</svg>
</span>
</div>


### Tetrahedral stereochemistry
Tetrahedral stereochemistry is denoted by a left parenthesis `(` meaning "counterclockwise"
or a right parenthesis `)` meaning "clockwise," referring to the first three neighbors of
a stereocenter atoms as they appear in the string, with the last neighbor (or implicit hydrogen)
in back.

| AMSR | molecule |
| --- | --- |
C(C.FO | [(S)-1-fluoroethanol](https://pubchem.ncbi.nlm.nih.gov/compound/57518764)
C)C.FO | [(1R)-1-fluoroethanol](https://pubchem.ncbi.nlm.nih.gov/compound/60205193)


<div>
<span style="margin:20px"><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='160px' height='160px' viewBox='0 0 160 160'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='160.0' height='160.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 56.5,102.1 L 57.2,103.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 46.2,107.2 L 47.6,109.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 36.0,112.4 L 38.0,115.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 25.7,117.5 L 28.5,122.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 15.5,122.6 L 18.9,128.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 5.2,127.7 L 9.3,134.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 66.8,97.0 L 66.8,68.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 66.8,68.3 L 66.8,39.7' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 66.8,97.0 L 91.2,111.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 91.2,111.1 L 115.6,125.2' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-2' d='M 61.0 18.6
L 72.5 18.6
L 72.5 20.8
L 63.6 20.8
L 63.6 26.7
L 71.6 26.7
L 71.6 28.9
L 63.6 28.9
L 63.6 38.0
L 61.0 38.0
L 61.0 18.6
' fill='#33CCCC'/>
<path  class='atom-3' d='M 117.3 131.4
Q 117.3 126.7, 119.6 124.1
Q 121.9 121.5, 126.2 121.5
Q 130.6 121.5, 132.9 124.1
Q 135.2 126.7, 135.2 131.4
Q 135.2 136.1, 132.8 138.8
Q 130.5 141.4, 126.2 141.4
Q 122.0 141.4, 119.6 138.8
Q 117.3 136.1, 117.3 131.4
M 126.2 139.2
Q 129.2 139.2, 130.8 137.3
Q 132.4 135.3, 132.4 131.4
Q 132.4 127.5, 130.8 125.6
Q 129.2 123.7, 126.2 123.7
Q 123.3 123.7, 121.7 125.6
Q 120.1 127.5, 120.1 131.4
Q 120.1 135.3, 121.7 137.3
Q 123.3 139.2, 126.2 139.2
' fill='#FF0000'/>
<path  class='atom-3' d='M 137.5 121.7
L 140.1 121.7
L 140.1 130.0
L 150.1 130.0
L 150.1 121.7
L 152.7 121.7
L 152.7 141.1
L 150.1 141.1
L 150.1 132.2
L 140.1 132.2
L 140.1 141.1
L 137.5 141.1
L 137.5 121.7
' fill='#FF0000'/>
</svg>
</span>
<span style="margin:20px"><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='160px' height='160px' viewBox='0 0 160 160'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='160.0' height='160.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 66.8,97.0 L 5.2,127.7 L 9.3,134.9 Z' style='fill:#000000;fill-rule:evenodd;fill-opacity:1;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='bond-1 atom-0 atom-2' d='M 66.8,97.0 L 66.8,68.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 66.8,68.3 L 66.8,39.7' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 66.8,97.0 L 91.2,111.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 91.2,111.1 L 115.6,125.2' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-2' d='M 61.0 18.6
L 72.5 18.6
L 72.5 20.8
L 63.6 20.8
L 63.6 26.7
L 71.6 26.7
L 71.6 28.9
L 63.6 28.9
L 63.6 38.0
L 61.0 38.0
L 61.0 18.6
' fill='#33CCCC'/>
<path  class='atom-3' d='M 117.3 131.4
Q 117.3 126.7, 119.6 124.1
Q 121.9 121.5, 126.2 121.5
Q 130.6 121.5, 132.9 124.1
Q 135.2 126.7, 135.2 131.4
Q 135.2 136.1, 132.8 138.8
Q 130.5 141.4, 126.2 141.4
Q 122.0 141.4, 119.6 138.8
Q 117.3 136.1, 117.3 131.4
M 126.2 139.2
Q 129.2 139.2, 130.8 137.3
Q 132.4 135.3, 132.4 131.4
Q 132.4 127.5, 130.8 125.6
Q 129.2 123.7, 126.2 123.7
Q 123.3 123.7, 121.7 125.6
Q 120.1 127.5, 120.1 131.4
Q 120.1 135.3, 121.7 137.3
Q 123.3 139.2, 126.2 139.2
' fill='#FF0000'/>
<path  class='atom-3' d='M 137.5 121.7
L 140.1 121.7
L 140.1 130.0
L 150.1 130.0
L 150.1 121.7
L 152.7 121.7
L 152.7 141.1
L 150.1 141.1
L 150.1 132.2
L 140.1 132.2
L 140.1 141.1
L 137.5 141.1
L 137.5 121.7
' fill='#FF0000'/>
</svg>
</span>
</div>


## Developing
This repo uses pre-commit, so after cloning run `pip install -r requirements.txt` and `pre-commit install` prior to committing.
