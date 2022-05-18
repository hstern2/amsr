# AMSR
**A**nother **M**olecular **S**tring **R**epresentation,
inspired by

- [H. Hiz, "A Linearization of Chemical Graphs," *J. Chem. Doc.* **4**, 173-180 (1964)](https://pubs.acs.org/doi/10.1021/c160014a015)
- [SMILES](https://pubs.acs.org/doi/10.1021/ci00057a005)
- [PATTY](https://pubs.acs.org/doi/10.1021/ci00015a015)
- [DeepSMILES](https://github.com/baoilleach/deepsmiles)
- [SELFIES](https://github.com/aspuru-guzik-group/selfies)

## Installing
```
pip install amsr
```

## Usage
```py
import amsr

amsr.ToMol("CNcncc5cNcN6C.oC.o") # caffeine
```
![caffine](https://user-images.githubusercontent.com/19351218/151638119-b1439d47-5e5a-417e-9254-c34568e2f3d1.png)

```py
taxol_smi = "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C"

amsr.FromSmiles(taxol_smi)
# CccCC(C(C(C)C(OC4.CC)C)6coC(8[OAc].C.O....[OAc].O[Bz].CC(6OcoC(O.C)N[Bz].[Ph]....O.C.C

amsr.FromSmiles(taxol_smi, useGroups=False)
# CccCC(C(C(C)C(OC4.CC)C)6coC(8OcoC..C.O....OcoC..Ococccccc6......CC(6OcoC(O.C)Ncocccccc6......cccccc6.........O.C.C
```

## Description

A molecular string representation in which every sequence of tokens generates a "reasonable" molecule.
You may have different ideas about what constitutes a reasonable molecule than this string representation.

### Atoms
Atoms are represented by their symbol
enclosed in square brackets, as in SMILES.
For a one-letter symbol,
brackets may be omitted.  Atoms are assumed to have a fixed
valence that limits the number of covalently-bonded neighbors.
If an atom makes fewer bonds than its valence, hydrogens are assumed.

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
explicitly like any other atom.  In the first
example below, the fluorines are added to the second carbon; the chlorine
is then added to the first carbon, since the second can no longer bond.
In the second example, after four fluorines are added, the carbon can
no longer bond, so a new molecule is formed.

| AMSR | molecule |
| --- | --- |
CCFFF[Cl] | [2-chloro-1,1,1-trifluoroethane](https://pubchem.ncbi.nlm.nih.gov/compound/2-Chloro-1_1_1-trifluoroethane)
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
<span style="margin:20px"><?xml version='1.0' encoding='iso-8859-1'?>
<svg version='1.1' baseProfile='full'
              xmlns='http://www.w3.org/2000/svg'
                      xmlns:rdkit='http://www.rdkit.org/xml'
                      xmlns:xlink='http://www.w3.org/1999/xlink'
                  xml:space='preserve'
width='160px' height='160px' viewBox='0 0 160 160'>
<!-- END OF HEADER -->
<rect style='opacity:1.0;fill:#FFFFFF;stroke:none' width='160.0' height='160.0' x='0.0' y='0.0'> </rect>
<path class='bond-0 atom-0 atom-1' d='M 80.0,80.0 L 108.4,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 108.4,80.0 L 136.7,80.0' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 80.0,80.0 L 51.6,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 51.6,80.0 L 23.3,80.0' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 80.0,80.0 L 80.0,53.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 80.0,53.5 L 80.0,26.9' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-0 atom-4' d='M 80.0,80.0 L 80.0,106.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-0 atom-4' d='M 80.0,106.5 L 80.0,133.1' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-1' d='M 138.3 71.0
L 149.1 71.0
L 149.1 73.0
L 140.8 73.0
L 140.8 78.5
L 148.2 78.5
L 148.2 80.6
L 140.8 80.6
L 140.8 89.0
L 138.3 89.0
L 138.3 71.0
' fill='#33CCCC'/>
<path  class='atom-2' d='M 10.9 71.0
L 21.7 71.0
L 21.7 73.0
L 13.4 73.0
L 13.4 78.5
L 20.7 78.5
L 20.7 80.6
L 13.4 80.6
L 13.4 89.0
L 10.9 89.0
L 10.9 71.0
' fill='#33CCCC'/>
<path  class='atom-3' d='M 74.6 7.3
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
<path  class='atom-4' d='M 74.6 134.7
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
</svg>
</span>
</div>


### Branches
Branches are formed automatically when atoms can no longer
make bonds.  They can also be made by "capping" or
"saturating" an atom with hydrogens, using a period `.`
(capping hydrogens are applied to the most
recently-added atom that can still make bonds).
New atoms will then be bonded to those added earlier, forming a branch.

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


### Rings
Rings are denoted by a single digit
(or two or more digits enclosed in square brackets)
giving the size of the ring.
A new bond is formed between
the two most recently-added atoms that
can make bonds and when bonded will form a ring of that size.

| AMSR | molecule |
| --- | --- |
CCO3 | [oxirane](https://en.wikipedia.org/wiki/Ethylene_oxide)
CCCCCC6 | [cyclohexane](https://en.wikipedia.org/wiki/Cyclohexane)
CCCCCCCCCCCC[12] | [cyclododecane](https://en.wikipedia.org/wiki/Cyclododecane)


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
</div>


### Double bonds (sp<sup>2</sup> centers)
Atoms making a double bond are indicated by changing
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


### Ring selection
When more than one ring of a given size can be formed, one or more single quote marks `'` immediately after
the digit will make ring-forming bonds with atoms appearing earlier in the
string, rather than the most recent.

| AMSR | molecule |
| --- | --- |
ccOcc5cccc6 | [benzofuran](https://en.wikipedia.org/wiki/Benzofuran)
ccOcc5cccc6' | [isobenzofuran](https://en.wikipedia.org/wiki/Isobenzofuran)


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
<path class='bond-0 atom-0 atom-1' d='M 33.4,116.0 L 7.3,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 36.7,105.4 L 18.4,80.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-0' d='M 75.7,102.2 L 33.4,116.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 7.3,80.0 L 17.6,65.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 17.6,65.8 L 27.8,51.7' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 40.3,46.3 L 58.0,52.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 58.0,52.0 L 75.7,57.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 75.7,57.8 L 75.7,102.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 66.8,64.4 L 66.8,95.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-8 atom-3' d='M 114.2,35.5 L 75.7,57.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-4 atom-5' d='M 75.7,102.2 L 114.2,124.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-5 atom-6' d='M 114.2,124.5 L 152.7,102.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-5 atom-6' d='M 115.5,113.4 L 142.5,97.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-6 atom-7' d='M 152.7,102.2 L 152.7,57.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-7 atom-8' d='M 152.7,57.8 L 114.2,35.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-7 atom-8' d='M 142.5,62.1 L 115.5,46.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-2' d='M 27.6 44.1
Q 27.6 41.0, 29.1 39.3
Q 30.6 37.7, 33.4 37.7
Q 36.2 37.7, 37.7 39.3
Q 39.2 41.0, 39.2 44.1
Q 39.2 47.1, 37.7 48.9
Q 36.2 50.6, 33.4 50.6
Q 30.6 50.6, 29.1 48.9
Q 27.6 47.1, 27.6 44.1
M 33.4 49.2
Q 35.3 49.2, 36.4 47.9
Q 37.4 46.6, 37.4 44.1
Q 37.4 41.6, 36.4 40.3
Q 35.3 39.1, 33.4 39.1
Q 31.5 39.1, 30.4 40.3
Q 29.4 41.6, 29.4 44.1
Q 29.4 46.6, 30.4 47.9
Q 31.5 49.2, 33.4 49.2
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
<path class='bond-0 atom-0 atom-1' d='M 78.6,58.6 L 38.0,45.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 69.9,64.8 L 41.4,55.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-0' d='M 78.6,101.4 L 78.6,58.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-8 atom-0' d='M 115.7,37.2 L 78.6,58.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 38.0,45.4 L 28.0,59.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 28.0,59.1 L 18.1,72.8' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 18.2,87.4 L 28.1,101.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 28.1,101.0 L 38.0,114.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 38.0,114.6 L 78.6,101.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 41.4,104.5 L 69.9,95.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-4 atom-5' d='M 78.6,101.4 L 115.7,122.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-5 atom-6' d='M 115.7,122.8 L 152.7,101.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-5 atom-6' d='M 117.0,112.2 L 142.9,97.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-6 atom-7' d='M 152.7,101.4 L 152.7,58.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-7 atom-8' d='M 152.7,58.6 L 115.7,37.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-7 atom-8' d='M 142.9,62.8 L 117.0,47.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-2' d='M 7.3 80.0
Q 7.3 77.1, 8.7 75.5
Q 10.1 73.9, 12.8 73.9
Q 15.5 73.9, 17.0 75.5
Q 18.4 77.1, 18.4 80.0
Q 18.4 83.0, 16.9 84.7
Q 15.5 86.3, 12.8 86.3
Q 10.2 86.3, 8.7 84.7
Q 7.3 83.0, 7.3 80.0
M 12.8 84.9
Q 14.7 84.9, 15.7 83.7
Q 16.7 82.5, 16.7 80.0
Q 16.7 77.7, 15.7 76.5
Q 14.7 75.2, 12.8 75.2
Q 11.0 75.2, 10.0 76.4
Q 9.0 77.6, 9.0 80.0
Q 9.0 82.5, 10.0 83.7
Q 11.0 84.9, 12.8 84.9
' fill='#FF0000'/>
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


### Formal charges, radical electrons, isotopes
Positive/negative formal charges are designated by one or more of `+`/`-`.
Radical electrons are denoted by one or more asterisks `*`.  An isotopic mass is denoted by a number prefix
before the atomic symbol (in which case square brackets must be used even for a one-letter symbol).

| AMSR | molecule |
| --- | --- |
[Mg++]S!!\:ooO-O- | [magnesium sulfate](https://en.wikipedia.org/wiki/Magnesium_sulfate)
CC.C.CCCCC.C.N6O* | [(2,2,6,6-Tetramethylpiperidin-1-yl)oxyl or TEMPO](https://en.wikipedia.org/wiki/TEMPO)
[2H]O[2H] | [heavy water](https://en.wikipedia.org/wiki/Heavy_water)


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
<path  class='atom-0' d='M 96.2 94.2
L 92.5 94.2
L 89.6 71.0
L 81.9 94.2
L 78.1 94.2
L 70.4 71.0
L 67.5 94.2
L 63.8 94.2
L 67.4 65.8
L 72.2 65.8
L 80.0 88.6
L 87.8 65.8
L 92.6 65.8
L 96.2 94.2
' fill='#000000'/>
<path  class='atom-0' d='M 116.7 73.6
L 116.7 95.4
Q 116.7 104.2, 107.3 104.2
Q 102.8 104.2, 98.8 102.1
L 100.1 99.4
Q 102.2 100.4, 103.8 100.8
Q 105.3 101.2, 107.3 101.2
Q 110.3 101.2, 111.6 99.9
Q 112.9 98.6, 112.9 95.4
L 112.9 91.8
Q 110.8 94.5, 107.1 94.5
Q 102.8 94.5, 100.4 91.8
Q 98.0 89.2, 98.0 84.2
Q 98.0 79.0, 100.8 76.2
Q 103.7 73.2, 108.6 73.2
Q 111.2 73.2, 113.6 74.0
L 113.7 73.6
L 116.7 73.6
M 107.3 91.5
Q 109.4 91.5, 110.9 90.4
Q 112.4 89.2, 112.9 87.2
L 112.9 76.8
Q 110.9 76.2, 108.6 76.2
Q 105.5 76.2, 103.7 78.3
Q 101.9 80.4, 101.9 84.2
Q 101.9 87.7, 103.3 89.6
Q 104.7 91.5, 107.3 91.5
' fill='#000000'/>
<path  class='atom-0' d='M 118.1 63.5
Q 118.8 61.8, 120.4 60.8
Q 122.0 59.8, 124.3 59.8
Q 127.1 59.8, 128.7 61.3
Q 130.3 62.8, 130.3 65.6
Q 130.3 68.3, 128.2 70.9
Q 126.2 73.5, 122.0 76.6
L 130.6 76.6
L 130.6 78.7
L 118.0 78.7
L 118.0 76.9
Q 121.5 74.4, 123.6 72.6
Q 125.6 70.7, 126.7 69.1
Q 127.7 67.4, 127.7 65.7
Q 127.7 63.9, 126.8 62.9
Q 125.9 61.9, 124.3 61.9
Q 122.8 61.9, 121.8 62.5
Q 120.8 63.1, 120.1 64.5
L 118.1 63.5
' fill='#000000'/>
<path  class='atom-0' d='M 133.0 69.6
L 138.0 69.6
L 138.0 64.4
L 140.2 64.4
L 140.2 69.6
L 145.3 69.6
L 145.3 71.5
L 140.2 71.5
L 140.2 76.8
L 138.0 76.8
L 138.0 71.5
L 133.0 71.5
L 133.0 69.6
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
<path class='bond-0 atom-0 atom-1' d='M 46.0,77.1 L 32.5,114.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 46.0,77.1 L 7.3,70.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 46.0,77.1 L 46.0,37.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-9 atom-0' d='M 74.8,93.7 L 60.4,85.4' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-9 atom-0' d='M 60.4,85.4 L 46.0,77.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 46.0,37.8 L 80.0,18.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 80.0,18.2 L 114.0,37.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-6' d='M 114.0,37.8 L 114.0,77.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-6 atom-7' d='M 114.0,77.1 L 152.7,70.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-6 atom-8' d='M 114.0,77.1 L 127.5,114.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-6 atom-9' d='M 114.0,77.1 L 99.6,85.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-6 atom-9' d='M 99.6,85.4 L 85.2,93.7' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 80.0,103.3 L 80.0,116.4' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 80.0,116.4 L 80.0,129.6' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-9' d='M 77.5 91.2
L 81.2 97.1
Q 81.5 97.7, 82.1 98.7
Q 82.7 99.8, 82.7 99.8
L 82.7 91.2
L 84.2 91.2
L 84.2 102.3
L 82.7 102.3
L 78.8 95.9
Q 78.3 95.1, 77.8 94.2
Q 77.4 93.4, 77.2 93.1
L 77.2 102.3
L 75.8 102.3
L 75.8 91.2
L 77.5 91.2
' fill='#0000FF'/>
<path  class='atom-10' d='M 74.9 136.1
Q 74.9 133.4, 76.2 131.9
Q 77.5 130.4, 80.0 130.4
Q 82.5 130.4, 83.8 131.9
Q 85.1 133.4, 85.1 136.1
Q 85.1 138.8, 83.8 140.3
Q 82.4 141.8, 80.0 141.8
Q 77.5 141.8, 76.2 140.3
Q 74.9 138.8, 74.9 136.1
M 80.0 140.6
Q 81.7 140.6, 82.6 139.5
Q 83.5 138.3, 83.5 136.1
Q 83.5 133.9, 82.6 132.8
Q 81.7 131.7, 80.0 131.7
Q 78.3 131.7, 77.4 132.8
Q 76.5 133.9, 76.5 136.1
Q 76.5 138.3, 77.4 139.5
Q 78.3 140.6, 80.0 140.6
' fill='#FF0000'/>
<path d='M 88.3,136.0 L 88.2,136.0 L 88.2,135.9 L 88.2,135.8 L 88.2,135.8 L 88.2,135.7 L 88.1,135.7 L 88.1,135.6 L 88.1,135.5 L 88.0,135.5 L 88.0,135.4 L 87.9,135.4 L 87.9,135.4 L 87.8,135.3 L 87.7,135.3 L 87.7,135.3 L 87.6,135.3 L 87.6,135.3 L 87.5,135.3 L 87.4,135.3 L 87.3,135.3 L 87.3,135.3 L 87.2,135.3 L 87.2,135.3 L 87.1,135.4 L 87.0,135.4 L 87.0,135.4 L 86.9,135.5 L 86.9,135.5 L 86.8,135.6 L 86.8,135.6 L 86.8,135.7 L 86.7,135.7 L 86.7,135.8 L 86.7,135.9 L 86.7,135.9 L 86.7,136.0 L 86.7,136.1 L 86.7,136.1 L 86.7,136.2 L 86.7,136.3 L 86.7,136.3 L 86.8,136.4 L 86.8,136.5 L 86.8,136.5 L 86.9,136.6 L 86.9,136.6 L 87.0,136.7 L 87.0,136.7 L 87.1,136.7 L 87.2,136.8 L 87.2,136.8 L 87.3,136.8 L 87.3,136.8 L 87.4,136.8 L 87.5,136.8 L 87.6,136.8 L 87.6,136.8 L 87.7,136.8 L 87.7,136.8 L 87.8,136.7 L 87.9,136.7 L 87.9,136.7 L 88.0,136.6 L 88.0,136.6 L 88.1,136.5 L 88.1,136.5 L 88.1,136.4 L 88.2,136.4 L 88.2,136.3 L 88.2,136.2 L 88.2,136.2 L 88.2,136.1 L 88.3,136.0 L 87.5,136.0 Z' style='fill:#000000;fill-rule:evenodd;fill-opacity:1;stroke:#000000;stroke-width:0.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
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
<path class='bond-0 atom-0 atom-1' d='M 135.8,69.6 L 115.8,81.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 115.8,81.1 L 95.9,92.6' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 74.8,92.5 L 54.9,81.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 54.9,81.0 L 34.9,69.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-0' d='M 126.5 53.8
Q 127.0 52.6, 128.1 51.9
Q 129.2 51.2, 130.8 51.2
Q 132.7 51.2, 133.8 52.3
Q 134.9 53.3, 134.9 55.2
Q 134.9 57.1, 133.5 58.9
Q 132.1 60.7, 129.2 62.8
L 135.1 62.8
L 135.1 64.2
L 126.5 64.2
L 126.5 63.0
Q 128.9 61.3, 130.3 60.0
Q 131.7 58.8, 132.4 57.6
Q 133.1 56.5, 133.1 55.3
Q 133.1 54.0, 132.5 53.4
Q 131.9 52.7, 130.8 52.7
Q 129.8 52.7, 129.1 53.1
Q 128.4 53.5, 127.9 54.4
L 126.5 53.8
' fill='#000000'/>
<path  class='atom-0' d='M 137.5 54.5
L 140.1 54.5
L 140.1 62.8
L 150.1 62.8
L 150.1 54.5
L 152.7 54.5
L 152.7 74.0
L 150.1 74.0
L 150.1 65.0
L 140.1 65.0
L 140.1 74.0
L 137.5 74.0
L 137.5 54.5
' fill='#000000'/>
<path  class='atom-1' d='M 76.6 98.7
Q 76.6 94.0, 78.9 91.4
Q 81.2 88.8, 85.5 88.8
Q 89.8 88.8, 92.1 91.4
Q 94.4 94.0, 94.4 98.7
Q 94.4 103.4, 92.1 106.1
Q 89.8 108.8, 85.5 108.8
Q 81.2 108.8, 78.9 106.1
Q 76.6 103.4, 76.6 98.7
M 85.5 106.6
Q 88.5 106.6, 90.1 104.6
Q 91.7 102.6, 91.7 98.7
Q 91.7 94.9, 90.1 92.9
Q 88.5 91.0, 85.5 91.0
Q 82.5 91.0, 80.9 92.9
Q 79.3 94.8, 79.3 98.7
Q 79.3 102.6, 80.9 104.6
Q 82.5 106.6, 85.5 106.6
' fill='#FF0000'/>
<path  class='atom-2' d='M 7.3 53.8
Q 7.8 52.6, 8.9 51.9
Q 10.0 51.2, 11.6 51.2
Q 13.5 51.2, 14.6 52.3
Q 15.7 53.3, 15.7 55.2
Q 15.7 57.1, 14.3 58.9
Q 12.9 60.7, 10.0 62.8
L 15.9 62.8
L 15.9 64.2
L 7.3 64.2
L 7.3 63.0
Q 9.7 61.3, 11.1 60.0
Q 12.5 58.8, 13.2 57.6
Q 13.9 56.5, 13.9 55.3
Q 13.9 54.0, 13.3 53.4
Q 12.7 52.7, 11.6 52.7
Q 10.6 52.7, 9.9 53.1
Q 9.2 53.5, 8.7 54.4
L 7.3 53.8
' fill='#000000'/>
<path  class='atom-2' d='M 18.3 54.5
L 20.9 54.5
L 20.9 62.8
L 30.9 62.8
L 30.9 54.5
L 33.5 54.5
L 33.5 74.0
L 30.9 74.0
L 30.9 65.0
L 20.9 65.0
L 20.9 74.0
L 18.3 74.0
L 18.3 54.5
' fill='#000000'/>
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
C(C.FO | [(1*S*)-1-fluoroethanol](https://pubchem.ncbi.nlm.nih.gov/compound/57518764)
C)C.FO | [(1*R*)-1-fluoroethanol](https://pubchem.ncbi.nlm.nih.gov/compound/60205193)


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


### *E/Z* stereochemistry
Stereochemistry for a double bond is denoted by an underscore `_` meaning "trans" or *E*,
or caret `^` meaning "cis" or *Z*, between the two atoms making the bond,
where the reference neighboring atoms are those that appear earliest in the string.

| AMSR | molecule |
| --- | --- |
c[Br][Cl]_c[Cl] | [(*E*)-1-bromo-1,2-dichloroethene](https://pubchem.ncbi.nlm.nih.gov/compound/E_-1-Bromo-1_2-dichloroethene)
c[Br][Cl]^c[Cl] | [(*Z*)-1-Bromo-1,2-dichloroethene](https://pubchem.ncbi.nlm.nih.gov/compound/Z_-1-Bromo-1_2-dichloroethene)


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
<path class='bond-0 atom-0 atom-1' d='M 61.9,92.1 L 43.4,102.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 43.4,102.8 L 24.8,113.5' style='fill:none;fill-rule:evenodd;stroke:#7F4C19;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 61.9,92.1 L 61.9,72.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 61.9,72.7 L 61.9,53.2' style='fill:none;fill-rule:evenodd;stroke:#00CC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 61.9,92.1 L 103.0,115.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 63.4,103.9 L 92.1,120.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 103.0,115.8 L 120.4,105.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 120.4,105.8 L 137.8,95.8' style='fill:none;fill-rule:evenodd;stroke:#00CC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-1' d='M 13.8 115.5
Q 15.1 115.9, 15.7 116.7
Q 16.4 117.4, 16.4 118.6
Q 16.4 120.5, 15.2 121.5
Q 14.0 122.6, 11.8 122.6
L 7.3 122.6
L 7.3 109.1
L 11.2 109.1
Q 13.5 109.1, 14.7 110.1
Q 15.9 111.0, 15.9 112.7
Q 15.9 114.7, 13.8 115.5
M 9.1 110.6
L 9.1 114.9
L 11.2 114.9
Q 12.6 114.9, 13.3 114.3
Q 14.0 113.8, 14.0 112.7
Q 14.0 110.6, 11.2 110.6
L 9.1 110.6
M 11.8 121.1
Q 13.1 121.1, 13.8 120.4
Q 14.5 119.8, 14.5 118.6
Q 14.5 117.5, 13.7 116.9
Q 13.0 116.4, 11.5 116.4
L 9.1 116.4
L 9.1 121.1
L 11.8 121.1
' fill='#7F4C19'/>
<path  class='atom-1' d='M 19.5 112.8
L 19.7 114.2
Q 20.7 112.6, 22.4 112.6
Q 22.9 112.6, 23.6 112.8
L 23.3 114.4
Q 22.5 114.2, 22.1 114.2
Q 21.3 114.2, 20.7 114.6
Q 20.2 114.9, 19.8 115.6
L 19.8 122.6
L 18.0 122.6
L 18.0 112.8
L 19.5 112.8
' fill='#7F4C19'/>
<path  class='atom-2' d='M 56.8 45.1
Q 56.8 41.8, 58.3 40.0
Q 59.9 38.2, 62.9 38.2
Q 65.6 38.2, 67.1 40.2
L 65.9 41.2
Q 64.8 39.8, 62.9 39.8
Q 60.8 39.8, 59.8 41.2
Q 58.7 42.5, 58.7 45.1
Q 58.7 47.8, 59.8 49.1
Q 60.9 50.5, 63.1 50.5
Q 64.6 50.5, 66.3 49.6
L 66.8 51.0
Q 66.1 51.5, 65.0 51.7
Q 64.0 52.0, 62.8 52.0
Q 59.9 52.0, 58.3 50.2
Q 56.8 48.4, 56.8 45.1
' fill='#00CC00'/>
<path  class='atom-2' d='M 68.8 37.4
L 70.5 37.4
L 70.5 51.8
L 68.8 51.8
L 68.8 37.4
' fill='#00CC00'/>
<path  class='atom-4' d='M 139.0 92.6
Q 139.0 89.2, 140.5 87.5
Q 142.1 85.7, 145.1 85.7
Q 147.9 85.7, 149.3 87.7
L 148.1 88.7
Q 147.0 87.3, 145.1 87.3
Q 143.1 87.3, 142.0 88.6
Q 140.9 90.0, 140.9 92.6
Q 140.9 95.2, 142.0 96.6
Q 143.1 98.0, 145.3 98.0
Q 146.8 98.0, 148.5 97.1
L 149.0 98.5
Q 148.3 99.0, 147.3 99.2
Q 146.2 99.5, 145.0 99.5
Q 142.1 99.5, 140.5 97.7
Q 139.0 95.9, 139.0 92.6
' fill='#00CC00'/>
<path  class='atom-4' d='M 151.0 84.9
L 152.7 84.9
L 152.7 99.3
L 151.0 99.3
L 151.0 84.9
' fill='#00CC00'/>
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
<path class='bond-0 atom-0 atom-1' d='M 102.8,91.6 L 102.8,71.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 102.8,71.7 L 102.8,51.8' style='fill:none;fill-rule:evenodd;stroke:#7F4C19;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 102.8,91.6 L 120.3,101.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 120.3,101.7 L 137.7,111.7' style='fill:none;fill-rule:evenodd;stroke:#00CC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 102.8,91.6 L 61.5,115.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 91.9,86.9 L 62.9,103.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 61.5,115.4 L 41.9,104.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 41.9,104.1 L 22.3,92.8' style='fill:none;fill-rule:evenodd;stroke:#00CC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-1' d='M 104.8 43.6
Q 106.1 43.9, 106.7 44.7
Q 107.4 45.5, 107.4 46.7
Q 107.4 48.5, 106.2 49.6
Q 105.0 50.7, 102.8 50.7
L 98.2 50.7
L 98.2 37.1
L 102.2 37.1
Q 104.5 37.1, 105.7 38.1
Q 106.9 39.0, 106.9 40.7
Q 106.9 42.8, 104.8 43.6
M 100.0 38.7
L 100.0 42.9
L 102.2 42.9
Q 103.6 42.9, 104.2 42.4
Q 104.9 41.8, 104.9 40.7
Q 104.9 38.7, 102.2 38.7
L 100.0 38.7
M 102.8 49.1
Q 104.1 49.1, 104.8 48.5
Q 105.5 47.9, 105.5 46.7
Q 105.5 45.6, 104.7 45.0
Q 104.0 44.4, 102.5 44.4
L 100.0 44.4
L 100.0 49.1
L 102.8 49.1
' fill='#7F4C19'/>
<path  class='atom-1' d='M 110.5 40.8
L 110.7 42.2
Q 111.7 40.7, 113.4 40.7
Q 113.9 40.7, 114.7 40.9
L 114.4 42.5
Q 113.5 42.3, 113.1 42.3
Q 112.3 42.3, 111.8 42.6
Q 111.2 42.9, 110.8 43.7
L 110.8 50.7
L 109.0 50.7
L 109.0 40.8
L 110.5 40.8
' fill='#7F4C19'/>
<path  class='atom-2' d='M 138.9 115.9
Q 138.9 112.6, 140.5 110.8
Q 142.1 109.0, 145.1 109.0
Q 147.8 109.0, 149.3 111.0
L 148.1 112.0
Q 147.0 110.6, 145.1 110.6
Q 143.0 110.6, 141.9 112.0
Q 140.9 113.3, 140.9 115.9
Q 140.9 118.6, 142.0 120.0
Q 143.1 121.3, 145.3 121.3
Q 146.8 121.3, 148.5 120.4
L 149.0 121.9
Q 148.3 122.3, 147.3 122.6
Q 146.2 122.9, 145.0 122.9
Q 142.1 122.9, 140.5 121.1
Q 138.9 119.3, 138.9 115.9
' fill='#00CC00'/>
<path  class='atom-2' d='M 151.0 108.2
L 152.7 108.2
L 152.7 122.7
L 151.0 122.7
L 151.0 108.2
' fill='#00CC00'/>
<path  class='atom-4' d='M 7.3 92.1
Q 7.3 88.7, 8.8 86.9
Q 10.4 85.2, 13.4 85.2
Q 16.2 85.2, 17.7 87.1
L 16.4 88.2
Q 15.3 86.7, 13.4 86.7
Q 11.4 86.7, 10.3 88.1
Q 9.2 89.5, 9.2 92.1
Q 9.2 94.7, 10.3 96.1
Q 11.5 97.5, 13.6 97.5
Q 15.1 97.5, 16.8 96.6
L 17.4 98.0
Q 16.7 98.5, 15.6 98.7
Q 14.5 99.0, 13.4 99.0
Q 10.4 99.0, 8.8 97.2
Q 7.3 95.4, 7.3 92.1
' fill='#00CC00'/>
<path  class='atom-4' d='M 19.3 84.4
L 21.1 84.4
L 21.1 98.8
L 19.3 98.8
L 19.3 84.4
' fill='#00CC00'/>
</svg>
</span>
</div>


### Groups
The following abbreviations may be used to represent various functional groups:
```py
[Ac], [Bn], [Boc], [Bz], [CCl3], [CF3], [CHO], [CN], [COO-], [COOEt], [COOH], [COOMe], [Cbz], [Et], [Ms], [NHAc], [NHMe], [NMe2], [NO2], [OAc], [OEt], [OMe], [Ph], [Piv], [Tf], [Tol], [Ts], [benzene], [iBu], [iPr], [nBu], [nPr], [sBu], [tBu]
```

| AMSR | molecule |
| --- | --- |
N+C)[Bn][COO-] | [L-phenylalanine](https://en.wikipedia.org/wiki/Phenylalanine)
cccccc6[iBu]..CC.[COOH] | [ibuprofen](https://en.wikipedia.org/wiki/Ibuprofen)
C[Ph][Et]coNcoNco6 | [phenobarbital](https://en.wikipedia.org/wiki/Phenobarbital)


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
<path class='bond-0 atom-1 atom-0' d='M 55.8,76.4 L 59.9,66.0 L 57.3,65.3 Z' style='fill:#000000;fill-rule:evenodd;fill-opacity:1;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='bond-0 atom-1 atom-0' d='M 59.9,66.0 L 58.8,54.2 L 64.1,55.5 Z' style='fill:#0000FF;fill-rule:evenodd;fill-opacity:1;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='bond-0 atom-1 atom-0' d='M 59.9,66.0 L 57.3,65.3 L 58.8,54.2 Z' style='fill:#0000FF;fill-rule:evenodd;fill-opacity:1;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='bond-1 atom-1 atom-2' d='M 55.8,76.4 L 74.9,95.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-1 atom-9' d='M 55.8,76.4 L 29.9,83.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 74.9,95.3 L 100.9,88.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 100.9,88.2 L 107.7,62.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 107.1,85.6 L 111.9,67.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-8 atom-3' d='M 120.0,107.1 L 100.9,88.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 107.7,62.1 L 133.6,55.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-6' d='M 133.6,55.0 L 152.7,73.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-6' d='M 132.7,61.7 L 146.1,74.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-6 atom-7' d='M 152.7,73.9 L 145.9,99.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-7 atom-8' d='M 145.9,99.9 L 120.0,107.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-7 atom-8' d='M 140.6,95.8 L 122.5,100.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 27.3,82.8 L 24.5,93.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 24.5,93.6 L 21.7,104.3' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 32.5,84.2 L 29.7,94.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 29.7,94.9 L 26.9,105.7' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-9 atom-11' d='M 29.9,83.5 L 22.4,76.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-9 atom-11' d='M 22.4,76.1 L 14.9,68.7' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-0' d='M 60.9 46.6
L 63.4 50.6
Q 63.7 51.0, 64.1 51.7
Q 64.5 52.4, 64.5 52.5
L 64.5 46.6
L 65.5 46.6
L 65.5 54.2
L 64.5 54.2
L 61.8 49.8
Q 61.5 49.3, 61.1 48.7
Q 60.8 48.1, 60.7 47.9
L 60.7 54.2
L 59.7 54.2
L 59.7 46.6
L 60.9 46.6
' fill='#0000FF'/>
<path  class='atom-0' d='M 66.4 46.6
L 67.5 46.6
L 67.5 49.8
L 71.3 49.8
L 71.3 46.6
L 72.4 46.6
L 72.4 54.2
L 71.3 54.2
L 71.3 50.7
L 67.5 50.7
L 67.5 54.2
L 66.4 54.2
L 66.4 46.6
' fill='#0000FF'/>
<path  class='atom-0' d='M 75.1 55.3
Q 75.6 55.5, 75.8 55.8
Q 76.1 56.1, 76.1 56.6
Q 76.1 57.0, 75.8 57.4
Q 75.6 57.7, 75.2 57.9
Q 74.8 58.1, 74.3 58.1
Q 73.8 58.1, 73.4 57.9
Q 73.0 57.7, 72.7 57.3
L 73.1 56.9
Q 73.4 57.2, 73.7 57.4
Q 73.9 57.5, 74.3 57.5
Q 74.8 57.5, 75.1 57.2
Q 75.3 57.0, 75.3 56.6
Q 75.3 56.1, 75.1 55.9
Q 74.8 55.6, 74.2 55.6
L 73.8 55.6
L 73.8 55.1
L 74.1 55.1
Q 74.7 55.1, 74.9 54.9
Q 75.2 54.6, 75.2 54.2
Q 75.2 53.9, 75.0 53.7
Q 74.8 53.5, 74.3 53.5
Q 73.9 53.5, 73.7 53.6
Q 73.4 53.8, 73.2 54.1
L 72.7 53.9
Q 72.9 53.5, 73.3 53.2
Q 73.7 52.9, 74.3 52.9
Q 75.1 52.9, 75.5 53.2
Q 75.9 53.6, 75.9 54.2
Q 75.9 54.6, 75.7 54.9
Q 75.5 55.2, 75.1 55.3
' fill='#0000FF'/>
<path  class='atom-0' d='M 73.0 47.9
L 74.4 47.9
L 74.4 46.5
L 75.0 46.5
L 75.0 47.9
L 76.3 47.9
L 76.3 48.4
L 75.0 48.4
L 75.0 49.9
L 74.4 49.9
L 74.4 48.4
L 73.0 48.4
L 73.0 47.9
' fill='#0000FF'/>
<path  class='atom-10' d='M 19.6 109.5
Q 19.6 107.7, 20.5 106.7
Q 21.4 105.7, 23.1 105.7
Q 24.8 105.7, 25.7 106.7
Q 26.6 107.7, 26.6 109.5
Q 26.6 111.4, 25.7 112.4
Q 24.8 113.5, 23.1 113.5
Q 21.4 113.5, 20.5 112.4
Q 19.6 111.4, 19.6 109.5
M 23.1 112.6
Q 24.3 112.6, 24.9 111.8
Q 25.5 111.1, 25.5 109.5
Q 25.5 108.0, 24.9 107.3
Q 24.3 106.5, 23.1 106.5
Q 21.9 106.5, 21.3 107.3
Q 20.7 108.0, 20.7 109.5
Q 20.7 111.1, 21.3 111.8
Q 21.9 112.6, 23.1 112.6
' fill='#FF0000'/>
<path  class='atom-11' d='M 7.3 64.6
Q 7.3 62.8, 8.2 61.8
Q 9.1 60.8, 10.8 60.8
Q 12.5 60.8, 13.4 61.8
Q 14.3 62.8, 14.3 64.6
Q 14.3 66.5, 13.3 67.5
Q 12.4 68.6, 10.8 68.6
Q 9.1 68.6, 8.2 67.5
Q 7.3 66.5, 7.3 64.6
M 10.8 67.7
Q 11.9 67.7, 12.6 66.9
Q 13.2 66.2, 13.2 64.6
Q 13.2 63.1, 12.6 62.4
Q 11.9 61.6, 10.8 61.6
Q 9.6 61.6, 9.0 62.4
Q 8.3 63.1, 8.3 64.6
Q 8.3 66.2, 9.0 66.9
Q 9.6 67.7, 10.8 67.7
' fill='#FF0000'/>
<path  class='atom-11' d='M 14.6 62.0
L 17.2 62.0
L 17.2 62.6
L 14.6 62.6
L 14.6 62.0
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
<path class='bond-0 atom-0 atom-1' d='M 91.5,85.8 L 71.1,87.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 88.2,81.9 L 73.9,82.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-0' d='M 100.6,67.5 L 91.5,85.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 71.1,87.0 L 59.8,70.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 59.8,70.0 L 69.0,51.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 64.9,69.1 L 71.3,56.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-2 atom-10' d='M 59.8,70.0 L 39.4,71.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 69.0,51.7 L 89.4,50.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 89.4,50.4 L 100.6,67.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 87.6,55.2 L 95.5,67.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-5 atom-6' d='M 100.6,67.5 L 121.0,66.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-6 atom-7' d='M 121.0,66.2 L 132.3,83.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-7 atom-8' d='M 132.3,83.3 L 152.7,82.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-7 atom-9' d='M 132.3,83.3 L 123.2,101.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-10 atom-11' d='M 39.4,71.2 L 28.2,54.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-12 atom-10 atom-12' d='M 39.4,71.2 L 30.3,89.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-13 atom-12 atom-13' d='M 30.2,87.5 L 21.6,88.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-13 atom-12 atom-13' d='M 21.6,88.0 L 13.0,88.5' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-13 atom-12 atom-13' d='M 30.5,91.6 L 21.8,92.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-13 atom-12 atom-13' d='M 21.8,92.1 L 13.2,92.6' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-14 atom-12 atom-14' d='M 30.3,89.5 L 34.8,96.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-14 atom-12 atom-14' d='M 34.8,96.3 L 39.3,103.1' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-13' d='M 7.3 90.8
Q 7.3 89.4, 8.0 88.6
Q 8.6 87.8, 9.9 87.8
Q 11.2 87.8, 11.9 88.6
Q 12.6 89.4, 12.6 90.8
Q 12.6 92.2, 11.9 93.0
Q 11.2 93.8, 9.9 93.8
Q 8.7 93.8, 8.0 93.0
Q 7.3 92.2, 7.3 90.8
M 9.9 93.1
Q 10.8 93.1, 11.3 92.5
Q 11.8 91.9, 11.8 90.8
Q 11.8 89.6, 11.3 89.1
Q 10.8 88.5, 9.9 88.5
Q 9.0 88.5, 8.6 89.1
Q 8.1 89.6, 8.1 90.8
Q 8.1 92.0, 8.6 92.5
Q 9.0 93.1, 9.9 93.1
' fill='#FF0000'/>
<path  class='atom-14' d='M 39.0 106.6
Q 39.0 105.2, 39.6 104.4
Q 40.3 103.6, 41.6 103.6
Q 42.9 103.6, 43.6 104.4
Q 44.3 105.2, 44.3 106.6
Q 44.3 108.0, 43.6 108.8
Q 42.9 109.6, 41.6 109.6
Q 40.3 109.6, 39.6 108.8
Q 39.0 108.0, 39.0 106.6
M 41.6 108.9
Q 42.5 108.9, 43.0 108.3
Q 43.5 107.7, 43.5 106.6
Q 43.5 105.4, 43.0 104.9
Q 42.5 104.3, 41.6 104.3
Q 40.7 104.3, 40.2 104.9
Q 39.8 105.4, 39.8 106.6
Q 39.8 107.7, 40.2 108.3
Q 40.7 108.9, 41.6 108.9
' fill='#FF0000'/>
<path  class='atom-14' d='M 45.0 103.7
L 45.8 103.7
L 45.8 106.2
L 48.7 106.2
L 48.7 103.7
L 49.5 103.7
L 49.5 109.5
L 48.7 109.5
L 48.7 106.8
L 45.8 106.8
L 45.8 109.5
L 45.0 109.5
L 45.0 103.7
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
<path class='bond-0 atom-0 atom-1' d='M 77.7,82.8 L 102.7,77.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-0 atom-7' d='M 77.7,82.8 L 87.0,106.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-0 atom-9' d='M 77.7,82.8 L 77.1,57.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-17 atom-15 atom-0' d='M 55.9,96.1 L 77.7,82.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 102.7,77.8 L 119.6,97.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 109.1,77.3 L 120.9,90.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-6 atom-1' d='M 110.9,53.6 L 102.7,77.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 119.6,97.0 L 144.6,91.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 144.6,91.9 L 152.7,67.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 141.0,86.7 L 146.7,69.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 152.7,67.8 L 135.9,48.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-6' d='M 135.9,48.6 L 110.9,53.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-6' d='M 133.1,54.4 L 115.6,57.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-7 atom-8' d='M 87.0,106.6 L 71.1,126.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 78.4,59.5 L 87.4,54.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 87.4,54.1 L 96.3,48.6' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 75.8,55.1 L 84.7,49.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 84.7,49.7 L 93.6,44.3' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-9 atom-11' d='M 77.1,57.3 L 67.6,52.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-9 atom-11' d='M 67.6,52.1 L 58.1,46.9' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-12 atom-11 atom-12' d='M 51.4,47.1 L 42.2,52.7' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-12 atom-11 atom-12' d='M 42.2,52.7 L 33.0,58.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-13 atom-12 atom-13' d='M 34.2,56.1 L 25.0,51.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-13 atom-12 atom-13' d='M 25.0,51.1 L 15.8,46.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-13 atom-12 atom-13' d='M 31.7,60.6 L 22.5,55.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-13 atom-12 atom-13' d='M 22.5,55.5 L 13.3,50.5' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-14 atom-12 atom-14' d='M 33.0,58.3 L 33.2,69.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-14 atom-12 atom-14' d='M 33.2,69.0 L 33.5,79.6' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-15 atom-14 atom-15' d='M 36.9,85.7 L 46.4,90.9' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-15 atom-14 atom-15' d='M 46.4,90.9 L 55.9,96.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-16 atom-15 atom-16' d='M 53.4,96.1 L 53.6,106.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-16 atom-15 atom-16' d='M 53.6,106.7 L 53.9,117.3' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-16 atom-15 atom-16' d='M 58.5,96.0 L 58.7,106.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-16 atom-15 atom-16' d='M 58.7,106.6 L 59.0,117.2' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-10' d='M 95.6 44.1
Q 95.6 42.3, 96.4 41.4
Q 97.3 40.4, 98.9 40.4
Q 100.5 40.4, 101.4 41.4
Q 102.2 42.3, 102.2 44.1
Q 102.2 45.8, 101.4 46.8
Q 100.5 47.8, 98.9 47.8
Q 97.3 47.8, 96.4 46.8
Q 95.6 45.8, 95.6 44.1
M 98.9 47.0
Q 100.0 47.0, 100.6 46.3
Q 101.2 45.5, 101.2 44.1
Q 101.2 42.7, 100.6 41.9
Q 100.0 41.2, 98.9 41.2
Q 97.8 41.2, 97.2 41.9
Q 96.6 42.6, 96.6 44.1
Q 96.6 45.5, 97.2 46.3
Q 97.8 47.0, 98.9 47.0
' fill='#FF0000'/>
<path  class='atom-11' d='M 53.2 41.5
L 55.5 45.3
Q 55.8 45.7, 56.1 46.4
Q 56.5 47.0, 56.5 47.1
L 56.5 41.5
L 57.5 41.5
L 57.5 48.7
L 56.5 48.7
L 54.0 44.5
Q 53.7 44.0, 53.3 43.5
Q 53.0 42.9, 52.9 42.7
L 52.9 48.7
L 52.0 48.7
L 52.0 41.5
L 53.2 41.5
' fill='#0000FF'/>
<path  class='atom-11' d='M 51.9 33.5
L 52.9 33.5
L 52.9 36.6
L 56.6 36.6
L 56.6 33.5
L 57.6 33.5
L 57.6 40.7
L 56.6 40.7
L 56.6 37.4
L 52.9 37.4
L 52.9 40.7
L 51.9 40.7
L 51.9 33.5
' fill='#0000FF'/>
<path  class='atom-13' d='M 7.3 46.1
Q 7.3 44.4, 8.1 43.4
Q 9.0 42.5, 10.6 42.5
Q 12.2 42.5, 13.0 43.4
Q 13.9 44.4, 13.9 46.1
Q 13.9 47.9, 13.0 48.9
Q 12.2 49.9, 10.6 49.9
Q 9.0 49.9, 8.1 48.9
Q 7.3 47.9, 7.3 46.1
M 10.6 49.1
Q 11.7 49.1, 12.3 48.3
Q 12.9 47.6, 12.9 46.1
Q 12.9 44.7, 12.3 44.0
Q 11.7 43.3, 10.6 43.3
Q 9.5 43.3, 8.9 44.0
Q 8.3 44.7, 8.3 46.1
Q 8.3 47.6, 8.9 48.3
Q 9.5 49.1, 10.6 49.1
' fill='#FF0000'/>
<path  class='atom-14' d='M 24.3 80.2
L 25.3 80.2
L 25.3 83.3
L 29.0 83.3
L 29.0 80.2
L 30.0 80.2
L 30.0 87.4
L 29.0 87.4
L 29.0 84.1
L 25.3 84.1
L 25.3 87.4
L 24.3 87.4
L 24.3 80.2
' fill='#0000FF'/>
<path  class='atom-14' d='M 32.0 80.2
L 34.3 84.0
Q 34.6 84.4, 34.9 85.1
Q 35.3 85.8, 35.3 85.8
L 35.3 80.2
L 36.3 80.2
L 36.3 87.4
L 35.3 87.4
L 32.8 83.3
Q 32.5 82.8, 32.2 82.2
Q 31.8 81.7, 31.8 81.5
L 31.8 87.4
L 30.8 87.4
L 30.8 80.2
L 32.0 80.2
' fill='#0000FF'/>
<path  class='atom-16' d='M 53.2 121.6
Q 53.2 119.9, 54.1 118.9
Q 54.9 117.9, 56.5 117.9
Q 58.1 117.9, 59.0 118.9
Q 59.8 119.9, 59.8 121.6
Q 59.8 123.3, 59.0 124.3
Q 58.1 125.3, 56.5 125.3
Q 54.9 125.3, 54.1 124.3
Q 53.2 123.3, 53.2 121.6
M 56.5 124.5
Q 57.6 124.5, 58.2 123.8
Q 58.8 123.0, 58.8 121.6
Q 58.8 120.2, 58.2 119.5
Q 57.6 118.7, 56.5 118.7
Q 55.4 118.7, 54.8 119.4
Q 54.2 120.2, 54.2 121.6
Q 54.2 123.0, 54.8 123.8
Q 55.4 124.5, 56.5 124.5
' fill='#FF0000'/>
</svg>
</span>
</div>


### Multiple molecules
More than one molecule may be specified by separating with `;`.

| AMSR | molecule |
| --- | --- |
ocC.O[benzene][COOH];CcoN[benzene]..O | [aspirin](https://en.wikipedia.org/wiki/Aspirin) and [acetaminophen](https://en.wikipedia.org/wiki/Paracetamol)


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
<path class='bond-0 atom-0 atom-1' d='M 31.9,97.7 L 33.8,106.5' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 33.8,106.5 L 35.7,115.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 36.2,96.8 L 38.0,105.6' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 38.0,105.6 L 39.9,114.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 37.8,114.8 L 21.6,129.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-1 atom-3' d='M 37.8,114.8 L 46.5,117.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-1 atom-3' d='M 46.5,117.6 L 55.1,120.4' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 61.9,118.5 L 68.3,112.7' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 68.3,112.7 L 74.7,106.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 74.7,106.9 L 70.1,85.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 78.2,102.8 L 75.0,87.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-9 atom-4' d='M 95.4,113.6 L 74.7,106.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-6' d='M 70.1,85.7 L 86.3,71.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-6 atom-7' d='M 86.3,71.1 L 107.0,77.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-6 atom-7' d='M 88.0,76.2 L 102.5,80.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-7 atom-8' d='M 107.0,77.8 L 111.5,99.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-8 atom-9' d='M 111.5,99.1 L 95.4,113.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-8 atom-9' d='M 106.2,98.0 L 94.9,108.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 95.4,113.6 L 99.9,134.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-10 atom-11' d='M 99.2,137.0 L 107.9,139.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-10 atom-11' d='M 107.9,139.8 L 116.6,142.6' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-10 atom-11' d='M 100.6,132.9 L 109.2,135.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-10 atom-11' d='M 109.2,135.7 L 117.9,138.5' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-12 atom-10 atom-12' d='M 99.9,134.9 L 93.5,140.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-12 atom-10 atom-12' d='M 93.5,140.7 L 87.1,146.5' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-13 atom-13 atom-14' d='M 15.5,37.0 L 36.6,31.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-14 atom-14 atom-15' d='M 38.7,32.0 L 41.0,23.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-14 atom-14 atom-15' d='M 41.0,23.3 L 43.3,14.7' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-14 atom-14 atom-15' d='M 34.5,30.9 L 36.8,22.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-14 atom-14 atom-15' d='M 36.8,22.2 L 39.1,13.6' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-15 atom-14 atom-16' d='M 36.6,31.4 L 42.8,37.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-15 atom-14 atom-16' d='M 42.8,37.7 L 49.1,43.9' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-16 atom-16 atom-17' d='M 54.8,46.0 L 63.9,43.6' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-16 atom-16 atom-17' d='M 63.9,43.6 L 73.0,41.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-17 atom-17 atom-18' d='M 73.0,41.2 L 78.6,20.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-17 atom-17 atom-18' d='M 78.0,39.1 L 82.0,24.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-22 atom-22 atom-17' d='M 88.4,56.6 L 73.0,41.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-18 atom-18 atom-19' d='M 78.6,20.1 L 99.6,14.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-19 atom-19 atom-20' d='M 99.6,14.5 L 115.0,29.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-19 atom-19 atom-20' d='M 98.9,19.9 L 109.6,30.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-20 atom-20 atom-21' d='M 115.0,29.9 L 109.4,50.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-23 atom-20 atom-23' d='M 115.0,29.9 L 123.9,27.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-23 atom-20 atom-23' d='M 123.9,27.5 L 132.7,25.2' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-21 atom-21 atom-22' d='M 109.4,50.9 L 88.4,56.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-21 atom-21 atom-22' d='M 105.1,47.6 L 90.4,51.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path  class='atom-0' d='M 30.4 93.6
Q 30.4 92.1, 31.1 91.2
Q 31.9 90.4, 33.2 90.4
Q 34.6 90.4, 35.3 91.2
Q 36.1 92.1, 36.1 93.6
Q 36.1 95.1, 35.3 95.9
Q 34.6 96.8, 33.2 96.8
Q 31.9 96.8, 31.1 95.9
Q 30.4 95.1, 30.4 93.6
M 33.2 96.1
Q 34.2 96.1, 34.7 95.4
Q 35.2 94.8, 35.2 93.6
Q 35.2 92.3, 34.7 91.7
Q 34.2 91.1, 33.2 91.1
Q 32.3 91.1, 31.8 91.7
Q 31.3 92.3, 31.3 93.6
Q 31.3 94.8, 31.8 95.4
Q 32.3 96.1, 33.2 96.1
' fill='#FF0000'/>
<path  class='atom-3' d='M 55.7 121.5
Q 55.7 120.1, 56.4 119.2
Q 57.1 118.4, 58.5 118.4
Q 59.9 118.4, 60.6 119.2
Q 61.3 120.1, 61.3 121.5
Q 61.3 123.0, 60.6 123.9
Q 59.8 124.7, 58.5 124.7
Q 57.1 124.7, 56.4 123.9
Q 55.7 123.0, 55.7 121.5
M 58.5 124.0
Q 59.4 124.0, 59.9 123.4
Q 60.5 122.8, 60.5 121.5
Q 60.5 120.3, 59.9 119.7
Q 59.4 119.1, 58.5 119.1
Q 57.6 119.1, 57.0 119.7
Q 56.5 120.3, 56.5 121.5
Q 56.5 122.8, 57.0 123.4
Q 57.6 124.0, 58.5 124.0
' fill='#FF0000'/>
<path  class='atom-11' d='M 117.8 141.7
Q 117.8 140.2, 118.5 139.3
Q 119.3 138.5, 120.6 138.5
Q 122.0 138.5, 122.7 139.3
Q 123.5 140.2, 123.5 141.7
Q 123.5 143.1, 122.7 144.0
Q 122.0 144.8, 120.6 144.8
Q 119.3 144.8, 118.5 144.0
Q 117.8 143.2, 117.8 141.7
M 120.6 144.2
Q 121.6 144.2, 122.1 143.5
Q 122.6 142.9, 122.6 141.7
Q 122.6 140.4, 122.1 139.8
Q 121.6 139.2, 120.6 139.2
Q 119.7 139.2, 119.2 139.8
Q 118.7 140.4, 118.7 141.7
Q 118.7 142.9, 119.2 143.5
Q 119.7 144.2, 120.6 144.2
' fill='#FF0000'/>
<path  class='atom-12' d='M 75.7 146.5
L 76.5 146.5
L 76.5 149.1
L 79.7 149.1
L 79.7 146.5
L 80.5 146.5
L 80.5 152.6
L 79.7 152.6
L 79.7 149.8
L 76.5 149.8
L 76.5 152.6
L 75.7 152.6
L 75.7 146.5
' fill='#FF0000'/>
<path  class='atom-12' d='M 80.9 149.5
Q 80.9 148.1, 81.7 147.2
Q 82.4 146.4, 83.8 146.4
Q 85.1 146.4, 85.9 147.2
Q 86.6 148.1, 86.6 149.5
Q 86.6 151.0, 85.8 151.9
Q 85.1 152.7, 83.8 152.7
Q 82.4 152.7, 81.7 151.9
Q 80.9 151.0, 80.9 149.5
M 83.8 152.0
Q 84.7 152.0, 85.2 151.4
Q 85.7 150.8, 85.7 149.5
Q 85.7 148.3, 85.2 147.7
Q 84.7 147.1, 83.8 147.1
Q 82.8 147.1, 82.3 147.7
Q 81.8 148.3, 81.8 149.5
Q 81.8 150.8, 82.3 151.4
Q 82.8 152.0, 83.8 152.0
' fill='#FF0000'/>
<path  class='atom-15' d='M 39.4 10.4
Q 39.4 8.9, 40.1 8.1
Q 40.8 7.3, 42.2 7.3
Q 43.6 7.3, 44.3 8.1
Q 45.0 8.9, 45.0 10.4
Q 45.0 11.9, 44.3 12.8
Q 43.6 13.6, 42.2 13.6
Q 40.8 13.6, 40.1 12.8
Q 39.4 11.9, 39.4 10.4
M 42.2 12.9
Q 43.1 12.9, 43.6 12.3
Q 44.2 11.6, 44.2 10.4
Q 44.2 9.2, 43.6 8.6
Q 43.1 8.0, 42.2 8.0
Q 41.3 8.0, 40.8 8.6
Q 40.2 9.2, 40.2 10.4
Q 40.2 11.7, 40.8 12.3
Q 41.3 12.9, 42.2 12.9
' fill='#FF0000'/>
<path  class='atom-16' d='M 50.6 43.7
L 52.6 47.0
Q 52.8 47.3, 53.1 47.9
Q 53.5 48.5, 53.5 48.5
L 53.5 43.7
L 54.3 43.7
L 54.3 49.9
L 53.5 49.9
L 51.3 46.3
Q 51.0 45.9, 50.8 45.4
Q 50.5 44.9, 50.4 44.8
L 50.4 49.9
L 49.6 49.9
L 49.6 43.7
L 50.6 43.7
' fill='#0000FF'/>
<path  class='atom-16' d='M 49.6 50.5
L 50.4 50.5
L 50.4 53.1
L 53.5 53.1
L 53.5 50.5
L 54.4 50.5
L 54.4 56.7
L 53.5 56.7
L 53.5 53.8
L 50.4 53.8
L 50.4 56.7
L 49.6 56.7
L 49.6 50.5
' fill='#0000FF'/>
<path  class='atom-23' d='M 133.2 24.3
Q 133.2 22.8, 134.0 22.0
Q 134.7 21.2, 136.1 21.2
Q 137.4 21.2, 138.2 22.0
Q 138.9 22.8, 138.9 24.3
Q 138.9 25.8, 138.2 26.6
Q 137.4 27.5, 136.1 27.5
Q 134.7 27.5, 134.0 26.6
Q 133.2 25.8, 133.2 24.3
M 136.1 26.8
Q 137.0 26.8, 137.5 26.2
Q 138.0 25.5, 138.0 24.3
Q 138.0 23.1, 137.5 22.5
Q 137.0 21.9, 136.1 21.9
Q 135.1 21.9, 134.6 22.5
Q 134.1 23.1, 134.1 24.3
Q 134.1 25.5, 134.6 26.2
Q 135.1 26.8, 136.1 26.8
' fill='#FF0000'/>
<path  class='atom-23' d='M 139.6 21.2
L 140.5 21.2
L 140.5 23.8
L 143.6 23.8
L 143.6 21.2
L 144.5 21.2
L 144.5 27.4
L 143.6 27.4
L 143.6 24.5
L 140.5 24.5
L 140.5 27.4
L 139.6 27.4
L 139.6 21.2
' fill='#FF0000'/>
</svg>
</span>
</div>


## Developing
This repo uses pre-commit, so after cloning run `pip install -r requirements.txt` and `pre-commit install` prior to committing.
