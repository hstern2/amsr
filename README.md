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
# CccCC`C`C`C'C`OC4.CC'C'6coC`8[OAc].C.O....[OAc].O[Bz].CC`6OcoC`O.C'N[Bz].[Ph]....O.C.C

amsr.FromSmiles(taxol_smi, useGroups=False)
# CccCC`C`C`C'C`OC4.CC'C'6coC`8OcoC..C.O....OcoC..Ococccccc6......CC`6OcoC`O.C'Ncocccccc6......cccccc6.........O.C.C
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
<path class='atom-0' d='M 46.1 73.0
Q 46.1 66.0, 49.4 62.3
Q 52.7 58.6, 59.0 58.6
Q 64.8 58.6, 68.0 62.7
L 65.3 64.9
Q 63.0 61.9, 59.0 61.9
Q 54.7 61.9, 52.4 64.8
Q 50.2 67.6, 50.2 73.0
Q 50.2 78.6, 52.5 81.5
Q 54.9 84.4, 59.4 84.4
Q 62.6 84.4, 66.2 82.5
L 67.3 85.5
Q 65.8 86.5, 63.6 87.0
Q 61.4 87.6, 58.9 87.6
Q 52.7 87.6, 49.4 83.8
Q 46.1 80.1, 46.1 73.0
' fill='#000000'/>
<path class='atom-0' d='M 71.7 58.9
L 75.6 58.9
L 75.6 71.0
L 90.1 71.0
L 90.1 58.9
L 93.9 58.9
L 93.9 87.2
L 90.1 87.2
L 90.1 74.2
L 75.6 74.2
L 75.6 87.2
L 71.7 87.2
L 71.7 58.9
' fill='#000000'/>
<path class='atom-0' d='M 111.6 95.0
L 113.9 95.0
L 113.9 97.1
L 111.6 97.1
L 111.6 101.4
L 109.2 101.4
L 109.2 97.1
L 99.5 97.1
L 99.5 95.4
L 107.7 82.7
L 111.6 82.7
L 111.6 95.0
M 102.6 95.0
L 109.2 95.0
L 109.2 84.4
L 102.6 95.0
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
<path class='atom-0' d='M 44.9 58.9
L 48.7 58.9
L 48.7 71.0
L 63.2 71.0
L 63.2 58.9
L 67.1 58.9
L 67.1 87.2
L 63.2 87.2
L 63.2 74.2
L 48.7 74.2
L 48.7 87.2
L 44.9 87.2
L 44.9 58.9
' fill='#FF0000'/>
<path class='atom-0' d='M 72.7 86.2
Q 73.4 84.5, 75.1 83.5
Q 76.7 82.5, 79.0 82.5
Q 81.8 82.5, 83.4 84.0
Q 85.0 85.6, 85.0 88.3
Q 85.0 91.1, 82.9 93.6
Q 80.9 96.2, 76.7 99.3
L 85.3 99.3
L 85.3 101.4
L 72.7 101.4
L 72.7 99.6
Q 76.2 97.1, 78.2 95.3
Q 80.3 93.5, 81.3 91.8
Q 82.3 90.1, 82.3 88.4
Q 82.3 86.6, 81.4 85.6
Q 80.5 84.6, 79.0 84.6
Q 77.5 84.6, 76.5 85.2
Q 75.5 85.8, 74.8 87.2
L 72.7 86.2
' fill='#FF0000'/>
<path class='atom-0' d='M 89.1 73.0
Q 89.1 66.2, 92.5 62.4
Q 95.8 58.6, 102.1 58.6
Q 108.4 58.6, 111.7 62.4
Q 115.1 66.2, 115.1 73.0
Q 115.1 79.9, 111.7 83.8
Q 108.3 87.7, 102.1 87.7
Q 95.9 87.7, 92.5 83.8
Q 89.1 79.9, 89.1 73.0
M 102.1 84.5
Q 106.4 84.5, 108.7 81.6
Q 111.1 78.7, 111.1 73.0
Q 111.1 67.4, 108.7 64.6
Q 106.4 61.8, 102.1 61.8
Q 97.8 61.8, 95.4 64.6
Q 93.1 67.4, 93.1 73.0
Q 93.1 78.7, 95.4 81.6
Q 97.8 84.5, 102.1 84.5
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
<path class='atom-0' d='M 52.1 66.7
L 55.9 66.7
L 55.9 78.7
L 70.4 78.7
L 70.4 66.7
L 74.2 66.7
L 74.2 95.0
L 70.4 95.0
L 70.4 81.9
L 55.9 81.9
L 55.9 95.0
L 52.1 95.0
L 52.1 66.7
' fill='#00CC00'/>
<path class='atom-0' d='M 79.2 80.8
Q 79.2 73.8, 82.5 70.1
Q 85.8 66.4, 92.1 66.4
Q 97.9 66.4, 101.1 70.5
L 98.4 72.6
Q 96.1 69.6, 92.1 69.6
Q 87.8 69.6, 85.5 72.5
Q 83.3 75.4, 83.3 80.8
Q 83.3 86.4, 85.6 89.3
Q 88.0 92.2, 92.5 92.2
Q 95.7 92.2, 99.3 90.3
L 100.4 93.3
Q 98.9 94.2, 96.7 94.8
Q 94.5 95.4, 92.0 95.4
Q 85.8 95.4, 82.5 91.6
Q 79.2 87.8, 79.2 80.8
' fill='#00CC00'/>
<path class='atom-0' d='M 104.2 64.6
L 107.9 64.6
L 107.9 95.0
L 104.2 95.0
L 104.2 64.6
' fill='#00CC00'/>
</svg>
</span>
</div>


### Chains
Each atom in a chain is bonded to the most recently added atom that
can still make bonds, according to its valence. Hydrogens may be added
explicitly like any other atom.  In the
example below, the fluorines are added to the second carbon; the chlorine
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
<path class='bond-0 atom-0 atom-1' d='M 63.1,104.7 L 105.8,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 105.8,80.0 L 124.3,69.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 124.3,69.3 L 142.8,58.7' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-1 atom-3' d='M 105.8,80.0 L 115.7,97.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-1 atom-3' d='M 115.7,97.1 L 125.5,114.1' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-1 atom-4' d='M 105.8,80.0 L 96.0,62.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-1 atom-4' d='M 96.0,62.9 L 86.1,45.9' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-0 atom-5' d='M 63.1,104.7 L 43.0,93.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-0 atom-5' d='M 43.0,93.1 L 22.9,81.5' style='fill:none;fill-rule:evenodd;stroke:#00CC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 65.2,103.4 L 63.1,104.7 L 62.1,104.1' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-2' d='M 144.4 48.3
L 152.7 48.3
L 152.7 49.9
L 146.3 49.9
L 146.3 54.2
L 152.0 54.2
L 152.0 55.8
L 146.3 55.8
L 146.3 62.3
L 144.4 62.3
L 144.4 48.3
' fill='#33CCCC'/>
<path class='atom-3' d='M 126.4 115.7
L 134.7 115.7
L 134.7 117.3
L 128.2 117.3
L 128.2 121.6
L 134.0 121.6
L 134.0 123.2
L 128.2 123.2
L 128.2 129.7
L 126.4 129.7
L 126.4 115.7
' fill='#33CCCC'/>
<path class='atom-4' d='M 77.0 30.3
L 85.3 30.3
L 85.3 31.9
L 78.9 31.9
L 78.9 36.1
L 84.6 36.1
L 84.6 37.7
L 78.9 37.7
L 78.9 44.3
L 77.0 44.3
L 77.0 30.3
' fill='#33CCCC'/>
<path class='atom-5' d='M 7.3 80.5
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
Q 7.3 84.0, 7.3 80.5
' fill='#00CC00'/>
<path class='atom-5' d='M 19.5 72.5
L 21.3 72.5
L 21.3 87.5
L 19.5 87.5
L 19.5 72.5
' fill='#00CC00'/>
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
<path d='M 53.3,53.4 L 55.8,52.0 L 58.2,53.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
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
<path class='bond-0 atom-0 atom-1' d='M 141.4,72.6 L 31.6,9.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 31.6,9.2 L 31.6,63.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 31.6,63.4 L 31.6,117.5' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-0' d='M 48.8,126.1 L 95.1,99.3' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-0' d='M 95.1,99.3 L 141.4,72.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 135.9,69.5 L 141.4,72.6 L 139.1,74.0' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 37.1,12.4 L 31.6,9.2 L 31.6,12.0' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-2' d='M 18.6 136.1
Q 18.6 129.3, 22.0 125.5
Q 25.3 121.7, 31.6 121.7
Q 37.9 121.7, 41.3 125.5
Q 44.6 129.3, 44.6 136.1
Q 44.6 143.0, 41.2 146.9
Q 37.8 150.8, 31.6 150.8
Q 25.4 150.8, 22.0 146.9
Q 18.6 143.0, 18.6 136.1
M 31.6 147.6
Q 35.9 147.6, 38.3 144.7
Q 40.6 141.8, 40.6 136.1
Q 40.6 130.5, 38.3 127.7
Q 35.9 124.9, 31.6 124.9
Q 27.3 124.9, 24.9 127.7
Q 22.6 130.5, 22.6 136.1
Q 22.6 141.8, 24.9 144.7
Q 27.3 147.6, 31.6 147.6
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
<path class='bond-1 atom-1 atom-2' d='M 116.4,143.0 L 43.6,143.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 43.6,143.0 L 7.3,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 7.3,80.0 L 43.6,17.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 43.6,17.0 L 116.4,17.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-0' d='M 116.4,17.0 L 152.7,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 150.9,83.1 L 152.7,80.0 L 150.9,76.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 118.2,139.8 L 116.4,143.0 L 112.7,143.0' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 47.3,143.0 L 43.6,143.0 L 41.8,139.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 9.1,83.1 L 7.3,80.0 L 9.1,76.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 41.8,20.2 L 43.6,17.0 L 47.3,17.0' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 112.7,17.0 L 116.4,17.0 L 118.2,20.2' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
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
<path class='bond-11 atom-11 atom-0' d='M 143.0,43.6 L 152.7,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 152.2,81.8 L 152.7,80.0 L 152.2,78.2' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 143.5,114.5 L 143.0,116.4 L 141.7,117.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 117.7,141.7 L 116.4,143.0 L 114.5,143.5' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 81.8,152.2 L 80.0,152.7 L 78.2,152.2' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 45.5,143.5 L 43.6,143.0 L 42.3,141.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 18.3,117.7 L 17.0,116.4 L 16.5,114.5' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 7.8,81.8 L 7.3,80.0 L 7.8,78.2' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 16.5,45.5 L 17.0,43.6 L 18.3,42.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 42.3,18.3 L 43.6,17.0 L 45.5,16.5' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 78.2,7.8 L 80.0,7.3 L 81.8,7.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 114.5,16.5 L 116.4,17.0 L 117.7,18.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 141.7,42.3 L 143.0,43.6 L 143.5,45.5' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
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
<path class='bond-0 atom-0 atom-1' d='M 9.1,70.1 L 64.9,70.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 64.9,70.1 L 120.6,70.1' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 9.1,89.4 L 64.9,89.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 64.9,89.4 L 120.6,89.4' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='atom-1' d='M 124.9 79.9
Q 124.9 73.1, 128.2 69.3
Q 131.6 65.5, 137.9 65.5
Q 144.1 65.5, 147.5 69.3
Q 150.9 73.1, 150.9 79.9
Q 150.9 86.7, 147.5 90.7
Q 144.1 94.5, 137.9 94.5
Q 131.6 94.5, 128.2 90.7
Q 124.9 86.8, 124.9 79.9
M 137.9 91.3
Q 142.2 91.3, 144.5 88.5
Q 146.9 85.5, 146.9 79.9
Q 146.9 74.3, 144.5 71.5
Q 142.2 68.7, 137.9 68.7
Q 133.5 68.7, 131.2 71.5
Q 128.9 74.3, 128.9 79.9
Q 128.9 85.6, 131.2 88.5
Q 133.5 91.3, 137.9 91.3
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
<path class='bond-0 atom-0 atom-1' d='M 140.1,80.0 L 110.1,132.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 116.4,143.0 L 43.6,143.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 43.6,143.0 L 7.3,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 49.9,132.1 L 19.9,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 7.3,80.0 L 43.6,17.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 43.6,17.0 L 116.4,17.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 49.9,27.9 L 110.1,27.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-0' d='M 116.4,17.0 L 152.7,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 150.9,83.1 L 152.7,80.0 L 150.9,76.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 118.2,139.8 L 116.4,143.0 L 112.7,143.0' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 47.3,143.0 L 43.6,143.0 L 41.8,139.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 9.1,83.1 L 7.3,80.0 L 9.1,76.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 41.8,20.2 L 43.6,17.0 L 47.3,17.0' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 112.7,17.0 L 116.4,17.0 L 118.2,20.2' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
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
<path class='bond-0 atom-0 atom-1' d='M 7.3,91.5 L 76.3,51.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 76.3,51.6 L 104.3,67.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 104.3,67.8 L 132.4,84.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 76.3,65.4 L 98.4,78.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 98.4,78.1 L 126.4,94.3' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 72.9,53.6 L 76.3,51.6 L 77.7,52.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-2' d='M 132.0 96.7
Q 132.0 91.3, 134.7 88.2
Q 137.4 85.2, 142.4 85.2
Q 147.4 85.2, 150.0 88.2
Q 152.7 91.3, 152.7 96.7
Q 152.7 102.2, 150.0 105.3
Q 147.3 108.4, 142.4 108.4
Q 137.4 108.4, 134.7 105.3
Q 132.0 102.2, 132.0 96.7
M 142.4 105.9
Q 145.8 105.9, 147.7 103.6
Q 149.5 101.2, 149.5 96.7
Q 149.5 92.3, 147.7 90.0
Q 145.8 87.8, 142.4 87.8
Q 138.9 87.8, 137.0 90.0
Q 135.2 92.2, 135.2 96.7
Q 135.2 101.3, 137.0 103.6
Q 138.9 105.9, 142.4 105.9
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
<path class='bond-0 atom-0 atom-1' d='M 128.1,73.9 L 90.2,21.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 95.0,7.3 L 16.6,32.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 16.6,32.7 L 16.6,115.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 29.0,41.7 L 29.0,106.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 16.6,115.1 L 49.1,125.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 49.1,125.7 L 81.5,136.2' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-0' d='M 105.5,126.1 L 124.4,100.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-0' d='M 124.4,100.0 L 143.4,73.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 141.0,70.6 L 143.4,73.9 L 142.4,75.2' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 97.4,10.6 L 95.0,7.3 L 91.0,8.5' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 20.5,31.5 L 16.6,32.7 L 16.6,36.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 16.6,111.0 L 16.6,115.1 L 18.2,115.6' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-4' d='M 84.3 140.6
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
<path class='bond-0 atom-0 atom-1' d='M 121.0,64.1 L 88.7,19.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 92.8,7.3 L 26.0,29.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 26.0,29.0 L 26.0,99.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 36.5,36.6 L 36.5,91.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 26.0,99.2 L 54.4,108.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 54.4,108.5 L 82.9,117.7' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-0' d='M 101.7,108.6 L 117.9,86.4' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-0' d='M 117.9,86.4 L 134.0,64.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 132.0,61.2 L 134.0,64.1 L 133.2,65.2' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 94.8,10.1 L 92.8,7.3 L 89.4,8.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 29.3,27.9 L 26.0,29.0 L 26.0,32.5' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 26.0,95.7 L 26.0,99.2 L 27.4,99.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-4' d='M 88.4 111.0
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
<path class='atom-4' d='M 85.0 132.8
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
When more than one ring of a given size can be formed, one or more `@` signs immediately after
the digit will make ring-forming bonds with atoms appearing earlier in the
string, rather than the most recent.

| AMSR | molecule |
| --- | --- |
ccOcc5cccc6 | [benzofuran](https://en.wikipedia.org/wiki/Benzofuran)
ccOcc5cccc6@ | [isobenzofuran](https://en.wikipedia.org/wiki/Isobenzofuran)


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
<path class='bond-0 atom-0 atom-1' d='M 36.0,108.1 L 15.5,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 7.3,80.0 L 17.4,66.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 17.4,66.0 L 27.6,52.1' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 40.7,46.4 L 58.2,52.1' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 58.2,52.1 L 75.7,57.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 75.7,57.8 L 75.7,102.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 69.0,62.6 L 69.0,97.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-0' d='M 75.7,102.2 L 33.4,116.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-4 atom-5' d='M 75.7,102.2 L 114.2,124.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-5 atom-6' d='M 114.2,124.5 L 152.7,102.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-5 atom-6' d='M 114.2,116.8 L 146.1,98.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-6 atom-7' d='M 152.7,102.2 L 152.7,57.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-7 atom-8' d='M 152.7,57.8 L 114.2,35.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-7 atom-8' d='M 146.1,61.6 L 114.2,43.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-8 atom-3' d='M 114.2,35.5 L 75.7,57.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 32.1,114.2 L 33.4,116.0 L 35.5,115.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 8.6,81.8 L 7.3,80.0 L 7.8,79.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 112.3,123.4 L 114.2,124.5 L 116.1,123.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 150.8,103.3 L 152.7,102.2 L 152.7,100.0' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 152.7,60.0 L 152.7,57.8 L 150.8,56.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 116.1,36.6 L 114.2,35.5 L 112.3,36.6' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-2' d='M 27.6 44.1
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
<path class='bond-0 atom-0 atom-1' d='M 72.2,63.3 L 40.4,52.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 38.0,45.4 L 28.1,58.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 28.1,58.9 L 18.3,72.5' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 18.4,87.7 L 28.2,101.2' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 28.2,101.2 L 38.0,114.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 38.0,114.6 L 78.6,101.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 40.4,107.1 L 72.2,96.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-0' d='M 78.6,101.4 L 78.6,58.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-4 atom-5' d='M 78.6,101.4 L 115.7,122.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-5 atom-6' d='M 115.7,122.8 L 152.7,101.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-5 atom-6' d='M 115.7,115.4 L 146.3,97.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-6 atom-7' d='M 152.7,101.4 L 152.7,58.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-7 atom-8' d='M 152.7,58.6 L 115.7,37.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-7 atom-8' d='M 146.3,62.3 L 115.7,44.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-8 atom-0' d='M 115.7,37.2 L 78.6,58.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 40.0,46.1 L 38.0,45.4 L 37.5,46.1' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 37.5,113.9 L 38.0,114.6 L 40.0,113.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 113.8,121.7 L 115.7,122.8 L 117.5,121.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 150.9,102.5 L 152.7,101.4 L 152.7,99.2' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 152.7,60.8 L 152.7,58.6 L 150.9,57.5' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 117.5,38.3 L 115.7,37.2 L 113.8,38.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-2' d='M 7.3 80.0
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
<path class='bond-0 atom-0 atom-1' d='M 151.0,80.0 L 92.9,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 92.9,80.0 L 34.8,80.0' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 151.0,99.7 L 92.9,99.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 92.9,99.7 L 34.8,99.7' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 151.0,60.3 L 92.9,60.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 92.9,60.3 L 34.8,60.3' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='atom-1' d='M 13.4 65.8
L 22.7 80.8
Q 23.6 82.3, 25.1 85.0
Q 26.6 87.7, 26.7 87.8
L 26.7 65.8
L 30.4 65.8
L 30.4 94.2
L 26.6 94.2
L 16.6 77.8
Q 15.4 75.8, 14.2 73.6
Q 13.0 71.4, 12.6 70.8
L 12.6 94.2
L 9.0 94.2
L 9.0 65.8
L 13.4 65.8
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
<path class='bond-0 atom-0 atom-1' d='M 133.9,84.7 L 111.5,84.7' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 111.5,84.7 L 89.2,84.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 133.9,75.0 L 111.5,75.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 111.5,75.0 L 89.2,75.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 70.8,84.7 L 48.5,84.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 48.5,84.7 L 26.1,84.7' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 70.8,75.0 L 48.5,75.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 48.5,75.0 L 26.1,75.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='atom-0' d='M 136.0 79.9
Q 136.0 75.5, 138.2 73.1
Q 140.3 70.6, 144.4 70.6
Q 148.4 70.6, 150.6 73.1
Q 152.7 75.5, 152.7 79.9
Q 152.7 84.3, 150.5 86.9
Q 148.4 89.4, 144.4 89.4
Q 140.3 89.4, 138.2 86.9
Q 136.0 84.4, 136.0 79.9
M 144.4 87.3
Q 147.1 87.3, 148.6 85.4
Q 150.2 83.6, 150.2 79.9
Q 150.2 76.3, 148.6 74.5
Q 147.1 72.7, 144.4 72.7
Q 141.6 72.7, 140.1 74.5
Q 138.6 76.3, 138.6 79.9
Q 138.6 83.6, 140.1 85.4
Q 141.6 87.3, 144.4 87.3
' fill='#FF0000'/>
<path class='atom-1' d='M 73.0 79.9
Q 73.0 75.4, 75.1 73.0
Q 77.2 70.6, 81.3 70.6
Q 85.0 70.6, 87.0 73.3
L 85.3 74.7
Q 83.9 72.8, 81.3 72.8
Q 78.5 72.8, 77.0 74.6
Q 75.6 76.4, 75.6 79.9
Q 75.6 83.5, 77.1 85.4
Q 78.6 87.2, 81.5 87.2
Q 83.6 87.2, 85.9 86.0
L 86.6 88.0
Q 85.7 88.6, 84.2 88.9
Q 82.8 89.3, 81.2 89.3
Q 77.2 89.3, 75.1 86.9
Q 73.0 84.5, 73.0 79.9
' fill='#000000'/>
<path class='atom-2' d='M 7.3 79.9
Q 7.3 75.5, 9.4 73.1
Q 11.6 70.6, 15.6 70.6
Q 19.7 70.6, 21.8 73.1
Q 24.0 75.5, 24.0 79.9
Q 24.0 84.3, 21.8 86.9
Q 19.6 89.4, 15.6 89.4
Q 11.6 89.4, 9.4 86.9
Q 7.3 84.4, 7.3 79.9
M 15.6 87.3
Q 18.4 87.3, 19.9 85.4
Q 21.4 83.6, 21.4 79.9
Q 21.4 76.3, 19.9 74.5
Q 18.4 72.7, 15.6 72.7
Q 12.9 72.7, 11.3 74.5
Q 9.8 76.3, 9.8 79.9
Q 9.8 83.6, 11.3 85.4
Q 12.9 87.3, 15.6 87.3
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
<path class='bond-0 atom-0 atom-1' d='M 7.3,106.9 L 38.4,88.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 38.4,88.9 L 69.6,70.9' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 90.4,70.9 L 121.6,88.9' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 121.6,88.9 L 152.7,106.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='atom-1' d='M 73.3 73.1
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
<path class='bond-0 atom-0 atom-1' d='M 7.3,149.0 L 38.4,131.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 38.4,131.0 L 69.6,113.0' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 86.3,92.4 L 86.3,65.3' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 86.3,65.3 L 86.3,38.2' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 73.7,92.4 L 73.7,65.3' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 73.7,65.3 L 73.7,38.2' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-1 atom-3' d='M 90.4,113.0 L 121.6,131.0' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-1 atom-3' d='M 121.6,131.0 L 152.7,149.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='atom-1' d='M 73.3 115.2
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
<path class='atom-2' d='M 69.1 23.1
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
<path class='bond-0 atom-0 atom-1' d='M 72.1,84.6 L 52.2,96.0' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 52.2,96.0 L 32.3,107.5' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 80.0,91.6 L 80.0,112.1' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 80.0,112.1 L 80.0,132.6' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 87.9,84.6 L 107.8,96.0' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 107.8,96.0 L 127.7,107.5' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-0 atom-4' d='M 80.0,68.9 L 80.0,48.2' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-0 atom-4' d='M 80.0,48.2 L 80.0,27.4' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-0 atom-5' d='M 87.9,75.4 L 107.8,64.0' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-0 atom-5' d='M 107.8,64.0 L 127.7,52.5' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-0 atom-6' d='M 72.1,75.4 L 52.2,64.0' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-0 atom-6' d='M 52.2,64.0 L 32.3,52.5' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='atom-0' d='M 74.9 86.2
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
<path class='atom-1' d='M 19.5 102.8
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
<path class='atom-2' d='M 74.6 134.7
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
<path class='atom-3' d='M 129.8 102.8
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
<path class='atom-4' d='M 74.6 7.3
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
<path class='atom-5' d='M 129.8 39.1
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
<path class='atom-6' d='M 19.5 39.1
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
<path class='atom-0' d='M 67.3 92.2
L 63.7 92.2
L 60.7 69.0
L 53.1 92.2
L 49.2 92.2
L 41.5 69.0
L 38.6 92.2
L 34.9 92.2
L 38.5 63.9
L 43.3 63.9
L 51.1 86.7
L 58.9 63.9
L 63.7 63.9
L 67.3 92.2
' fill='#000000'/>
<path class='atom-0' d='M 90.8 71.6
L 90.8 93.4
Q 90.8 102.2, 81.4 102.2
Q 76.9 102.2, 72.9 100.1
L 74.2 97.4
Q 76.3 98.5, 77.9 98.9
Q 79.4 99.3, 81.4 99.3
Q 84.4 99.3, 85.7 97.9
Q 87.0 96.6, 87.0 93.5
L 87.0 89.9
Q 84.9 92.5, 81.2 92.5
Q 76.9 92.5, 74.5 89.9
Q 72.1 87.2, 72.1 82.2
Q 72.1 77.1, 74.9 74.2
Q 77.8 71.3, 82.7 71.3
Q 85.3 71.3, 87.7 72.0
L 87.8 71.6
L 90.8 71.6
M 81.4 89.5
Q 83.5 89.5, 85.0 88.4
Q 86.5 87.3, 87.0 85.2
L 87.0 74.8
Q 85.0 74.2, 82.7 74.2
Q 79.6 74.2, 77.8 76.3
Q 76.0 78.4, 76.0 82.2
Q 76.0 85.7, 77.4 87.7
Q 78.8 89.5, 81.4 89.5
' fill='#000000'/>
<path class='atom-0' d='M 96.4 61.5
Q 97.1 59.8, 98.7 58.8
Q 100.4 57.8, 102.6 57.8
Q 105.4 57.8, 107.0 59.3
Q 108.6 60.9, 108.6 63.6
Q 108.6 66.3, 106.6 68.9
Q 104.5 71.5, 100.3 74.6
L 108.9 74.6
L 108.9 76.7
L 96.3 76.7
L 96.3 74.9
Q 99.8 72.4, 101.9 70.6
Q 104.0 68.7, 105.0 67.1
Q 106.0 65.4, 106.0 63.7
Q 106.0 61.9, 105.1 60.9
Q 104.2 59.9, 102.6 59.9
Q 101.1 59.9, 100.1 60.5
Q 99.1 61.1, 98.4 62.5
L 96.4 61.5
' fill='#000000'/>
<path class='atom-0' d='M 112.7 67.6
L 117.7 67.6
L 117.7 62.4
L 119.9 62.4
L 119.9 67.6
L 125.1 67.6
L 125.1 69.5
L 119.9 69.5
L 119.9 74.8
L 117.7 74.8
L 117.7 69.5
L 112.7 69.5
L 112.7 67.6
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
<path class='bond-3 atom-3 atom-4' d='M 46.0,37.8 L 80.0,18.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 80.0,18.2 L 114.0,37.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-6' d='M 114.0,37.8 L 114.0,77.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-6 atom-7' d='M 114.0,77.1 L 152.7,70.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-6 atom-8' d='M 114.0,77.1 L 127.5,114.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-6 atom-9' d='M 114.0,77.1 L 99.8,85.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-6 atom-9' d='M 99.8,85.3 L 85.5,93.6' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-9 atom-0' d='M 74.5,93.6 L 60.2,85.3' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-9 atom-0' d='M 60.2,85.3 L 46.0,77.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 80.0,103.6 L 80.0,116.4' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 80.0,116.4 L 80.0,129.1' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 46.0,39.8 L 46.0,37.8 L 47.7,36.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 78.3,19.1 L 80.0,18.2 L 81.7,19.1' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 112.3,36.8 L 114.0,37.8 L 114.0,39.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-9' d='M 77.5 91.2
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
<path class='atom-10' d='M 74.9 136.1
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
<path class='atom-10' d='M 88.3,136.0 L 88.2,136.1 L 88.2,136.2 L 88.2,136.2 L 88.2,136.3 L 88.2,136.4 L 88.1,136.4 L 88.1,136.5 L 88.1,136.5 L 88.0,136.6 L 88.0,136.6 L 87.9,136.7 L 87.9,136.7 L 87.8,136.7 L 87.7,136.8 L 87.7,136.8 L 87.6,136.8 L 87.6,136.8 L 87.5,136.8 L 87.4,136.8 L 87.3,136.8 L 87.3,136.8 L 87.2,136.8 L 87.2,136.8 L 87.1,136.7 L 87.0,136.7 L 87.0,136.7 L 86.9,136.6 L 86.9,136.6 L 86.8,136.5 L 86.8,136.5 L 86.8,136.4 L 86.7,136.3 L 86.7,136.3 L 86.7,136.2 L 86.7,136.1 L 86.7,136.1 L 86.7,136.0 L 86.7,135.9 L 86.7,135.9 L 86.7,135.8 L 86.7,135.7 L 86.8,135.7 L 86.8,135.6 L 86.8,135.6 L 86.9,135.5 L 86.9,135.5 L 87.0,135.4 L 87.0,135.4 L 87.1,135.4 L 87.2,135.3 L 87.2,135.3 L 87.3,135.3 L 87.3,135.3 L 87.4,135.3 L 87.5,135.3 L 87.6,135.3 L 87.6,135.3 L 87.7,135.3 L 87.7,135.3 L 87.8,135.3 L 87.9,135.4 L 87.9,135.4 L 88.0,135.4 L 88.0,135.5 L 88.1,135.5 L 88.1,135.6 L 88.1,135.7 L 88.2,135.7 L 88.2,135.8 L 88.2,135.8 L 88.2,135.9 L 88.2,136.0 L 88.3,136.0 L 87.5,136.0 Z' style='fill:#000000;fill-rule:evenodd;fill-opacity:1;stroke:#000000;stroke-width:0.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
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
<path class='bond-0 atom-0 atom-1' d='M 112.8,82.0 L 87.4,82.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 87.4,82.0 L 62.1,82.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='atom-0' d='M 116.0 67.5
Q 116.7 65.8, 118.3 64.9
Q 119.8 63.9, 122.0 63.9
Q 124.7 63.9, 126.2 65.4
Q 127.7 66.8, 127.7 69.4
Q 127.7 72.1, 125.8 74.6
Q 123.8 77.0, 119.8 80.0
L 128.0 80.0
L 128.0 82.0
L 116.0 82.0
L 116.0 80.3
Q 119.3 77.9, 121.3 76.1
Q 123.3 74.4, 124.3 72.8
Q 125.2 71.2, 125.2 69.6
Q 125.2 67.8, 124.4 66.9
Q 123.5 65.9, 122.0 65.9
Q 120.6 65.9, 119.6 66.5
Q 118.6 67.1, 118.0 68.4
L 116.0 67.5
' fill='#000000'/>
<path class='atom-0' d='M 131.5 68.4
L 135.2 68.4
L 135.2 80.0
L 149.1 80.0
L 149.1 68.4
L 152.7 68.4
L 152.7 95.5
L 149.1 95.5
L 149.1 83.0
L 135.2 83.0
L 135.2 95.5
L 131.5 95.5
L 131.5 68.4
' fill='#000000'/>
<path class='atom-1' d='M 7.3 68.6
L 10.9 68.6
L 10.9 80.1
L 24.8 80.1
L 24.8 68.6
L 28.5 68.6
L 28.5 95.7
L 24.8 95.7
L 24.8 83.2
L 10.9 83.2
L 10.9 95.7
L 7.3 95.7
L 7.3 68.6
' fill='#FF0000'/>
<path class='atom-1' d='M 34.0 82.1
Q 34.0 75.6, 37.2 71.9
Q 40.5 68.3, 46.5 68.3
Q 52.5 68.3, 55.7 71.9
Q 58.9 75.6, 58.9 82.1
Q 58.9 88.6, 55.6 92.4
Q 52.4 96.1, 46.5 96.1
Q 40.5 96.1, 37.2 92.4
Q 34.0 88.7, 34.0 82.1
M 46.5 93.0
Q 50.6 93.0, 52.8 90.3
Q 55.1 87.5, 55.1 82.1
Q 55.1 76.7, 52.8 74.1
Q 50.6 71.3, 46.5 71.3
Q 42.3 71.3, 40.1 74.0
Q 37.9 76.7, 37.9 82.1
Q 37.9 87.5, 40.1 90.3
Q 42.3 93.0, 46.5 93.0
' fill='#FF0000'/>
</svg>
</span>
</div>


### Tetrahedral stereochemistry
Tetrahedral stereochemistry is denoted by a single quote `'` meaning "clockwise"
or a backtick `` ` `` meaning "counterclockwise," referring to the first three neighbors of
a stereocenter atoms as they appear in the string, with the last neighbor (or implicit hydrogen)
in back.

| AMSR | molecule |
| --- | --- |
C`C.FO | [(1*S*)-1-fluoroethanol](https://pubchem.ncbi.nlm.nih.gov/compound/57518764)
C'C.FO | [(1*R*)-1-fluoroethanol](https://pubchem.ncbi.nlm.nih.gov/compound/60205193)


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
<path class='bond-0 atom-0 atom-1' d='M 106.1,66.0 L 69.0,94.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 106.1,66.0 L 124.7,74.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 124.7,74.0 L 143.4,81.9' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-1 atom-3' d='M 69.0,94.0 L 51.4,86.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-1 atom-3' d='M 51.4,86.5 L 33.9,79.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 104.2,67.4 L 106.1,66.0 L 107.0,66.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 70.9,92.6 L 69.0,94.0 L 68.1,93.6' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-2' d='M 144.9 77.7
L 152.7 77.7
L 152.7 79.2
L 146.7 79.2
L 146.7 83.1
L 152.1 83.1
L 152.1 84.7
L 146.7 84.7
L 146.7 90.8
L 144.9 90.8
L 144.9 77.7
' fill='#33CCCC'/>
<path class='atom-3' d='M 7.3 69.3
L 9.1 69.3
L 9.1 74.9
L 15.8 74.9
L 15.8 69.3
L 17.6 69.3
L 17.6 82.4
L 15.8 82.4
L 15.8 76.4
L 9.1 76.4
L 9.1 82.4
L 7.3 82.4
L 7.3 69.3
' fill='#FF0000'/>
<path class='atom-3' d='M 20.3 75.8
Q 20.3 72.7, 21.8 70.9
Q 23.4 69.1, 26.3 69.1
Q 29.2 69.1, 30.8 70.9
Q 32.3 72.7, 32.3 75.8
Q 32.3 79.0, 30.7 80.8
Q 29.2 82.6, 26.3 82.6
Q 23.4 82.6, 21.8 80.8
Q 20.3 79.0, 20.3 75.8
M 26.3 81.1
Q 28.3 81.1, 29.4 79.8
Q 30.5 78.4, 30.5 75.8
Q 30.5 73.2, 29.4 71.9
Q 28.3 70.6, 26.3 70.6
Q 24.3 70.6, 23.2 71.9
Q 22.1 73.2, 22.1 75.8
Q 22.1 78.5, 23.2 79.8
Q 24.3 81.1, 26.3 81.1
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
<path class='bond-0 atom-0 atom-1' d='M 106.1,66.0 L 69.0,94.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 106.1,66.0 L 124.7,74.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 124.7,74.0 L 143.4,81.9' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-1 atom-3' d='M 69.0,94.0 L 51.4,86.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-1 atom-3' d='M 51.4,86.5 L 33.9,79.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 104.2,67.4 L 106.1,66.0 L 107.0,66.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 70.9,92.6 L 69.0,94.0 L 68.1,93.6' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-2' d='M 144.9 77.7
L 152.7 77.7
L 152.7 79.2
L 146.7 79.2
L 146.7 83.1
L 152.1 83.1
L 152.1 84.7
L 146.7 84.7
L 146.7 90.8
L 144.9 90.8
L 144.9 77.7
' fill='#33CCCC'/>
<path class='atom-3' d='M 7.3 69.3
L 9.1 69.3
L 9.1 74.9
L 15.8 74.9
L 15.8 69.3
L 17.6 69.3
L 17.6 82.4
L 15.8 82.4
L 15.8 76.4
L 9.1 76.4
L 9.1 82.4
L 7.3 82.4
L 7.3 69.3
' fill='#FF0000'/>
<path class='atom-3' d='M 20.3 75.8
Q 20.3 72.7, 21.8 70.9
Q 23.4 69.1, 26.3 69.1
Q 29.2 69.1, 30.8 70.9
Q 32.3 72.7, 32.3 75.8
Q 32.3 79.0, 30.7 80.8
Q 29.2 82.6, 26.3 82.6
Q 23.4 82.6, 21.8 80.8
Q 20.3 79.0, 20.3 75.8
M 26.3 81.1
Q 28.3 81.1, 29.4 79.8
Q 30.5 78.4, 30.5 75.8
Q 30.5 73.2, 29.4 71.9
Q 28.3 70.6, 26.3 70.6
Q 24.3 70.6, 23.2 71.9
Q 22.1 73.2, 22.1 75.8
Q 22.1 78.5, 23.2 79.8
Q 24.3 81.1, 26.3 81.1
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
<path class='bond-0 atom-0 atom-1' d='M 62.0,92.1 L 43.6,102.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 43.6,102.8 L 25.2,113.4' style='fill:none;fill-rule:evenodd;stroke:#7F4C19;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 62.0,92.1 L 62.0,72.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 62.0,72.8 L 62.0,53.5' style='fill:none;fill-rule:evenodd;stroke:#00CC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 62.0,92.1 L 103.2,115.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 69.2,88.0 L 103.2,107.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 103.2,115.9 L 120.4,106.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 120.4,106.0 L 137.7,96.0' style='fill:none;fill-rule:evenodd;stroke:#00CC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 101.2,114.7 L 103.2,115.9 L 104.1,115.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-1' d='M 13.8 115.6
Q 15.1 115.9, 15.8 116.7
Q 16.4 117.5, 16.4 118.7
Q 16.4 120.5, 15.2 121.6
Q 14.0 122.6, 11.8 122.6
L 7.3 122.6
L 7.3 109.2
L 11.2 109.2
Q 13.6 109.2, 14.7 110.1
Q 15.9 111.0, 15.9 112.8
Q 15.9 114.8, 13.8 115.6
M 9.1 110.7
L 9.1 114.9
L 11.2 114.9
Q 12.6 114.9, 13.3 114.4
Q 14.0 113.8, 14.0 112.8
Q 14.0 110.7, 11.2 110.7
L 9.1 110.7
M 11.8 121.1
Q 13.1 121.1, 13.8 120.5
Q 14.5 119.9, 14.5 118.7
Q 14.5 117.6, 13.7 117.0
Q 13.0 116.4, 11.5 116.4
L 9.1 116.4
L 9.1 121.1
L 11.8 121.1
' fill='#7F4C19'/>
<path class='atom-1' d='M 19.5 112.9
L 19.7 114.2
Q 20.7 112.7, 22.4 112.7
Q 22.9 112.7, 23.6 112.9
L 23.4 114.5
Q 22.5 114.3, 22.1 114.3
Q 21.3 114.3, 20.8 114.6
Q 20.2 114.9, 19.8 115.7
L 19.8 122.6
L 18.0 122.6
L 18.0 112.9
L 19.5 112.9
' fill='#7F4C19'/>
<path class='atom-2' d='M 56.8 45.0
Q 56.8 41.7, 58.4 39.9
Q 60.0 38.2, 63.0 38.2
Q 65.7 38.2, 67.2 40.1
L 66.0 41.2
Q 64.9 39.7, 63.0 39.7
Q 60.9 39.7, 59.8 41.1
Q 58.8 42.4, 58.8 45.0
Q 58.8 47.7, 59.9 49.1
Q 61.0 50.4, 63.2 50.4
Q 64.7 50.4, 66.4 49.5
L 66.9 51.0
Q 66.2 51.4, 65.2 51.7
Q 64.1 52.0, 62.9 52.0
Q 60.0 52.0, 58.4 50.2
Q 56.8 48.4, 56.8 45.0
' fill='#00CC00'/>
<path class='atom-2' d='M 68.6 37.4
L 70.3 37.4
L 70.3 51.8
L 68.6 51.8
L 68.6 37.4
' fill='#00CC00'/>
<path class='atom-4' d='M 139.2 92.6
Q 139.2 89.3, 140.8 87.5
Q 142.4 85.7, 145.3 85.7
Q 148.1 85.7, 149.6 87.7
L 148.4 88.7
Q 147.3 87.3, 145.3 87.3
Q 143.3 87.3, 142.2 88.7
Q 141.2 90.0, 141.2 92.6
Q 141.2 95.3, 142.3 96.6
Q 143.4 98.0, 145.6 98.0
Q 147.0 98.0, 148.8 97.1
L 149.3 98.5
Q 148.6 99.0, 147.5 99.3
Q 146.5 99.5, 145.3 99.5
Q 142.4 99.5, 140.8 97.7
Q 139.2 95.9, 139.2 92.6
' fill='#00CC00'/>
<path class='atom-4' d='M 151.0 84.9
L 152.7 84.9
L 152.7 99.4
L 151.0 99.4
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
<path class='bond-0 atom-0 atom-1' d='M 102.9,91.6 L 102.9,71.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 102.9,71.9 L 102.9,52.1' style='fill:none;fill-rule:evenodd;stroke:#7F4C19;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 102.9,91.6 L 120.2,101.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 120.2,101.6 L 137.6,111.6' style='fill:none;fill-rule:evenodd;stroke:#00CC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 102.9,91.6 L 61.4,115.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 95.7,87.5 L 61.4,107.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 61.4,115.6 L 41.9,104.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 41.9,104.3 L 22.4,93.1' style='fill:none;fill-rule:evenodd;stroke:#00CC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 63.5,114.4 L 61.4,115.6 L 60.5,115.0' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-1' d='M 104.9 43.4
Q 106.2 43.8, 106.8 44.6
Q 107.5 45.4, 107.5 46.5
Q 107.5 48.4, 106.3 49.5
Q 105.1 50.5, 102.9 50.5
L 98.3 50.5
L 98.3 37.0
L 102.3 37.0
Q 104.6 37.0, 105.8 37.9
Q 107.0 38.9, 107.0 40.6
Q 107.0 42.6, 104.9 43.4
M 100.1 38.5
L 100.1 42.8
L 102.3 42.8
Q 103.6 42.8, 104.3 42.2
Q 105.0 41.7, 105.0 40.6
Q 105.0 38.5, 102.3 38.5
L 100.1 38.5
M 102.9 49.0
Q 104.2 49.0, 104.9 48.4
Q 105.6 47.7, 105.6 46.5
Q 105.6 45.4, 104.8 44.9
Q 104.0 44.3, 102.6 44.3
L 100.1 44.3
L 100.1 49.0
L 102.9 49.0
' fill='#7F4C19'/>
<path class='atom-1' d='M 110.6 40.7
L 110.8 42.1
Q 111.8 40.5, 113.5 40.5
Q 114.1 40.5, 114.8 40.7
L 114.5 42.3
Q 113.7 42.1, 113.2 42.1
Q 112.4 42.1, 111.9 42.5
Q 111.4 42.8, 110.9 43.5
L 110.9 50.5
L 109.1 50.5
L 109.1 40.7
L 110.6 40.7
' fill='#7F4C19'/>
<path class='atom-2' d='M 139.1 116.0
Q 139.1 112.7, 140.7 110.9
Q 142.3 109.1, 145.3 109.1
Q 148.1 109.1, 149.6 111.1
L 148.3 112.1
Q 147.2 110.7, 145.3 110.7
Q 143.3 110.7, 142.2 112.1
Q 141.1 113.4, 141.1 116.0
Q 141.1 118.7, 142.2 120.1
Q 143.3 121.5, 145.5 121.5
Q 147.0 121.5, 148.7 120.6
L 149.3 122.0
Q 148.6 122.5, 147.5 122.7
Q 146.4 123.0, 145.2 123.0
Q 142.3 123.0, 140.7 121.2
Q 139.1 119.4, 139.1 116.0
' fill='#00CC00'/>
<path class='atom-2' d='M 151.0 108.3
L 152.7 108.3
L 152.7 122.8
L 151.0 122.8
L 151.0 108.3
' fill='#00CC00'/>
<path class='atom-4' d='M 7.3 92.1
Q 7.3 88.7, 8.8 87.0
Q 10.4 85.2, 13.4 85.2
Q 16.2 85.2, 17.7 87.2
L 16.5 88.2
Q 15.4 86.8, 13.4 86.8
Q 11.4 86.8, 10.3 88.1
Q 9.2 89.5, 9.2 92.1
Q 9.2 94.8, 10.3 96.2
Q 11.5 97.5, 13.6 97.5
Q 15.1 97.5, 16.9 96.6
L 17.4 98.1
Q 16.7 98.5, 15.6 98.8
Q 14.6 99.1, 13.4 99.1
Q 10.4 99.1, 8.8 97.3
Q 7.3 95.5, 7.3 92.1
' fill='#00CC00'/>
<path class='atom-4' d='M 19.1 84.4
L 20.9 84.4
L 20.9 98.9
L 19.1 98.9
L 19.1 84.4
' fill='#00CC00'/>
</svg>
</span>
</div>


### Groups
The following abbreviations may be used to represent various functional groups:
```py
[@], [Ac], [Bn], [Boc], [Bz], [CCl3], [CF3], [CHO], [CN], [COO-], [COOEt], [COOH], [COOMe], [Cbz], [Cy], [Et], [Ms], [NC], [NHAc], [NHMe], [NMe2], [NO2], [OAc], [OEt], [OMe], [OiBu], [PO3], [Ph], [Pip], [Piv], [SMe], [SO3], [Tf], [Tol], [Ts], [iBu], [iPr], [nBu], [nDec], [nHept], [nHex], [nNon], [nOct], [nPent], [nPr], [sBu], [tBu]
```

| AMSR | molecule |
| --- | --- |
N+C'[Bn][COO-] | [L-phenylalanine](https://en.wikipedia.org/wiki/Phenylalanine)
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
<path class='bond-0 atom-1 atom-0' d='M 55.8,76.4 L 57.5,65.9 L 59.5,66.4 Z' style='fill:#000000;fill-rule:evenodd;fill-opacity:1;stroke:#000000;stroke-width:0.5px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='bond-0 atom-1 atom-0' d='M 57.5,65.9 L 63.1,56.5 L 59.2,55.4 Z' style='fill:#0000FF;fill-rule:evenodd;fill-opacity:1;stroke:#0000FF;stroke-width:0.5px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='bond-0 atom-1 atom-0' d='M 57.5,65.9 L 59.5,66.4 L 63.1,56.5 Z' style='fill:#0000FF;fill-rule:evenodd;fill-opacity:1;stroke:#0000FF;stroke-width:0.5px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='bond-1 atom-1 atom-2' d='M 55.8,76.4 L 74.9,95.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 74.9,95.3 L 100.9,88.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 100.9,88.2 L 107.7,62.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 105.4,86.9 L 111.0,65.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 107.7,62.1 L 133.6,55.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-6' d='M 133.6,55.0 L 152.7,73.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-6' d='M 132.4,59.5 L 148.2,75.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-6 atom-7' d='M 152.7,73.9 L 145.9,99.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-7 atom-8' d='M 145.9,99.9 L 120.0,107.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-7 atom-8' d='M 142.6,96.7 L 121.2,102.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-8 atom-3' d='M 120.0,107.1 L 100.9,88.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-1 atom-9' d='M 55.8,76.4 L 29.9,83.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 28.2,81.9 L 25.3,93.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 25.3,93.1 L 22.4,104.3' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 32.1,82.9 L 29.2,94.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 29.2,94.1 L 26.3,105.3' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-9 atom-11' d='M 29.9,83.5 L 22.5,76.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-9 atom-11' d='M 22.5,76.2 L 15.2,68.9' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 56.8,77.3 L 55.8,76.4 L 54.5,76.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 74.0,94.3 L 74.9,95.3 L 76.2,94.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 107.3,63.4 L 107.7,62.1 L 109.0,61.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 132.3,55.4 L 133.6,55.0 L 134.6,56.0' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 151.8,73.0 L 152.7,73.9 L 152.4,75.2' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 146.3,98.6 L 145.9,99.9 L 144.6,100.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 121.3,106.7 L 120.0,107.1 L 119.0,106.1' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 31.2,83.1 L 29.9,83.5 L 29.5,83.1' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-0' d='M 60.9 46.6
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
<path class='atom-0' d='M 66.9 46.6
L 68.0 46.6
L 68.0 49.8
L 71.9 49.8
L 71.9 46.6
L 72.9 46.6
L 72.9 54.2
L 71.9 54.2
L 71.9 50.7
L 68.0 50.7
L 68.0 54.2
L 66.9 54.2
L 66.9 46.6
' fill='#0000FF'/>
<path class='atom-0' d='M 76.8 55.3
Q 77.3 55.5, 77.5 55.8
Q 77.7 56.1, 77.7 56.6
Q 77.7 57.0, 77.5 57.4
Q 77.3 57.7, 76.9 57.9
Q 76.5 58.1, 76.0 58.1
Q 75.5 58.1, 75.1 57.9
Q 74.7 57.7, 74.4 57.3
L 74.8 56.9
Q 75.1 57.2, 75.3 57.4
Q 75.6 57.5, 76.0 57.5
Q 76.5 57.5, 76.8 57.2
Q 77.0 57.0, 77.0 56.6
Q 77.0 56.1, 76.7 55.9
Q 76.5 55.6, 75.8 55.6
L 75.5 55.6
L 75.5 55.1
L 75.8 55.1
Q 76.3 55.1, 76.6 54.9
Q 76.9 54.6, 76.9 54.2
Q 76.9 53.9, 76.7 53.7
Q 76.4 53.5, 76.0 53.5
Q 75.6 53.5, 75.4 53.6
Q 75.1 53.8, 74.9 54.1
L 74.4 53.9
Q 74.6 53.5, 75.0 53.2
Q 75.4 52.9, 76.0 52.9
Q 76.8 52.9, 77.2 53.2
Q 77.6 53.6, 77.6 54.2
Q 77.6 54.6, 77.4 54.9
Q 77.2 55.2, 76.8 55.3
' fill='#0000FF'/>
<path class='atom-0' d='M 74.4 47.9
L 75.7 47.9
L 75.7 46.5
L 76.3 46.5
L 76.3 47.9
L 77.7 47.9
L 77.7 48.4
L 76.3 48.4
L 76.3 49.9
L 75.7 49.9
L 75.7 48.4
L 74.4 48.4
L 74.4 47.9
' fill='#0000FF'/>
<path class='atom-10' d='M 19.6 109.5
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
<path class='atom-11' d='M 7.3 64.6
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
<path class='atom-11' d='M 15.3 62.0
L 18.0 62.0
L 18.0 62.6
L 15.3 62.6
L 15.3 62.0
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
<path class='bond-0 atom-0 atom-1' d='M 89.6,82.8 L 72.7,83.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 71.1,87.0 L 59.8,70.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 59.8,70.0 L 69.0,51.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 63.4,69.8 L 70.9,54.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 69.0,51.7 L 89.4,50.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 89.4,50.4 L 100.6,67.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 87.8,53.6 L 97.1,67.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-0' d='M 100.6,67.5 L 91.5,85.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-5 atom-6' d='M 100.6,67.5 L 121.0,66.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-6 atom-7' d='M 121.0,66.2 L 132.3,83.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-7 atom-8' d='M 132.3,83.3 L 152.7,82.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-7 atom-9' d='M 132.3,83.3 L 123.2,101.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-2 atom-10' d='M 59.8,70.0 L 39.4,71.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-10 atom-11' d='M 39.4,71.2 L 28.2,54.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-12 atom-10 atom-12' d='M 39.4,71.2 L 30.3,89.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-13 atom-12 atom-13' d='M 31.1,87.9 L 22.1,88.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-13 atom-12 atom-13' d='M 22.1,88.5 L 13.2,89.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-13 atom-12 atom-13' d='M 31.3,91.0 L 22.3,91.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-13 atom-12 atom-13' d='M 22.3,91.5 L 13.4,92.1' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-14 atom-12 atom-14' d='M 30.3,89.5 L 34.8,96.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-14 atom-12 atom-14' d='M 34.8,96.2 L 39.2,103.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 90.5,85.8 L 91.5,85.8 L 92.0,84.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 72.1,86.9 L 71.1,87.0 L 70.6,86.2' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 68.5,52.6 L 69.0,51.7 L 70.0,51.6' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 88.3,50.5 L 89.4,50.4 L 89.9,51.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 120.0,66.3 L 121.0,66.2 L 121.6,67.1' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 30.8,88.6 L 30.3,89.5 L 30.6,89.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-13' d='M 7.3 90.8
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
<path class='atom-14' d='M 39.0 106.6
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
<path class='atom-14' d='M 45.2 103.7
L 46.0 103.7
L 46.0 106.2
L 48.9 106.2
L 48.9 103.7
L 49.7 103.7
L 49.7 109.5
L 48.9 109.5
L 48.9 106.8
L 46.0 106.8
L 46.0 109.5
L 45.2 109.5
L 45.2 103.7
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
<path class='bond-1 atom-1 atom-2' d='M 102.7,77.8 L 119.6,97.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 107.0,76.9 L 121.0,92.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 119.6,97.0 L 144.6,91.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 144.6,91.9 L 152.7,67.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 141.6,88.6 L 148.4,68.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 152.7,67.8 L 135.9,48.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-6' d='M 135.9,48.6 L 110.9,53.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-6' d='M 134.5,52.8 L 113.8,57.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-6 atom-1' d='M 110.9,53.6 L 102.7,77.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-0 atom-7' d='M 77.7,82.8 L 87.0,106.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-7 atom-8' d='M 87.0,106.6 L 71.1,126.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-0 atom-9' d='M 77.7,82.8 L 77.1,57.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 77.2,59.5 L 86.5,53.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 86.5,53.9 L 95.7,48.2' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 75.2,56.3 L 84.5,50.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 84.5,50.6 L 93.8,45.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-9 atom-11' d='M 77.1,57.3 L 67.7,52.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-9 atom-11' d='M 67.7,52.2 L 58.3,47.0' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-12 atom-11 atom-12' d='M 51.2,47.3 L 42.1,52.8' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-12 atom-11 atom-12' d='M 42.1,52.8 L 33.0,58.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-13 atom-12 atom-13' d='M 34.8,57.2 L 25.3,51.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-13 atom-12 atom-13' d='M 25.3,51.9 L 15.7,46.7' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-13 atom-12 atom-13' d='M 33.0,60.5 L 23.4,55.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-13 atom-12 atom-13' d='M 23.4,55.3 L 13.8,50.1' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-14 atom-12 atom-14' d='M 33.0,58.3 L 33.2,68.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-14 atom-12 atom-14' d='M 33.2,68.9 L 33.5,79.4' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-15 atom-14 atom-15' d='M 37.1,85.8 L 46.5,90.9' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-15 atom-14 atom-15' d='M 46.5,90.9 L 55.9,96.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-16 atom-15 atom-16' d='M 54.0,95.0 L 54.2,106.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-16 atom-15 atom-16' d='M 54.2,106.1 L 54.5,117.1' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-16 atom-15 atom-16' d='M 57.8,94.9 L 58.1,106.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-16 atom-15 atom-16' d='M 58.1,106.0 L 58.3,117.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-17 atom-15 atom-0' d='M 55.9,96.1 L 77.7,82.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 118.7,96.0 L 119.6,97.0 L 120.8,96.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 143.3,92.2 L 144.6,91.9 L 145.0,90.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 152.3,69.0 L 152.7,67.8 L 151.9,66.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 136.7,49.6 L 135.9,48.6 L 134.6,48.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 112.1,53.4 L 110.9,53.6 L 110.5,54.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 86.5,105.4 L 87.0,106.6 L 86.2,107.6' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 77.2,58.6 L 77.1,57.3 L 76.7,57.1' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 33.4,58.1 L 33.0,58.3 L 33.0,58.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 55.5,95.8 L 55.9,96.1 L 57.0,95.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-10' d='M 95.6 44.1
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
<path class='atom-11' d='M 53.2 41.5
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
<path class='atom-11' d='M 51.9 33.5
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
<path class='atom-13' d='M 7.3 46.1
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
<path class='atom-14' d='M 23.7 80.2
L 24.7 80.2
L 24.7 83.3
L 28.4 83.3
L 28.4 80.2
L 29.4 80.2
L 29.4 87.4
L 28.4 87.4
L 28.4 84.1
L 24.7 84.1
L 24.7 87.4
L 23.7 87.4
L 23.7 80.2
' fill='#0000FF'/>
<path class='atom-14' d='M 32.0 80.2
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
<path class='atom-16' d='M 53.2 121.6
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
ocC.O[@][COOH];CcoN[@]..O | [aspirin](https://en.wikipedia.org/wiki/Aspirin) and [acetaminophen](https://en.wikipedia.org/wiki/Paracetamol)


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
<path class='bond-0 atom-0 atom-1' d='M 35.6,97.1 L 37.5,106.3' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 37.5,106.3 L 39.5,115.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 32.4,97.8 L 34.3,106.9' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 34.3,106.9 L 36.3,116.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 37.7,114.8 L 21.5,129.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-1 atom-3' d='M 37.7,114.8 L 46.3,117.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-1 atom-3' d='M 46.3,117.6 L 54.8,120.4' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 61.9,118.3 L 68.2,112.6' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 68.2,112.6 L 74.5,106.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 74.5,106.9 L 70.0,85.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 77.3,104.4 L 73.6,86.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-6' d='M 70.0,85.7 L 86.2,71.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-6 atom-7' d='M 86.2,71.1 L 106.9,77.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-6 atom-7' d='M 86.9,74.8 L 104.1,80.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-7 atom-8' d='M 106.9,77.8 L 111.4,99.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-8 atom-9' d='M 111.4,99.1 L 95.3,113.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-8 atom-9' d='M 107.8,97.9 L 94.5,110.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-9 atom-4' d='M 95.3,113.6 L 74.5,106.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 95.3,113.6 L 99.8,134.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-10 atom-11' d='M 98.4,136.2 L 107.4,139.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-10 atom-11' d='M 107.4,139.1 L 116.5,142.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-10 atom-11' d='M 99.4,133.1 L 108.4,136.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-10 atom-11' d='M 108.4,136.0 L 117.5,138.9' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-12 atom-10 atom-12' d='M 99.8,134.9 L 93.5,140.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-12 atom-10 atom-12' d='M 93.5,140.6 L 87.2,146.3' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-13 atom-13 atom-14' d='M 15.4,37.0 L 36.5,31.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-14 atom-14 atom-15' d='M 37.8,32.7 L 40.2,23.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-14 atom-14 atom-15' d='M 40.2,23.7 L 42.6,14.7' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-14 atom-14 atom-15' d='M 34.6,31.9 L 37.1,22.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-14 atom-14 atom-15' d='M 37.1,22.9 L 39.5,13.9' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-15 atom-14 atom-16' d='M 36.5,31.4 L 42.6,37.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-15 atom-14 atom-16' d='M 42.6,37.6 L 48.8,43.7' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-16 atom-16 atom-17' d='M 54.9,46.0 L 63.9,43.6' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-16 atom-16 atom-17' d='M 63.9,43.6 L 72.9,41.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-17 atom-17 atom-18' d='M 72.9,41.2 L 78.5,20.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-17 atom-17 atom-18' d='M 76.5,40.2 L 81.2,22.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-18 atom-18 atom-19' d='M 78.5,20.1 L 99.5,14.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-19 atom-19 atom-20' d='M 99.5,14.5 L 114.9,29.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-19 atom-19 atom-20' d='M 98.6,18.2 L 111.3,30.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-20 atom-20 atom-21' d='M 114.9,29.9 L 109.3,50.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-21 atom-21 atom-22' d='M 109.3,50.9 L 88.3,56.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-21 atom-21 atom-22' d='M 106.6,48.3 L 89.2,52.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-22 atom-22 atom-17' d='M 88.3,56.6 L 72.9,41.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-23 atom-20 atom-23' d='M 114.9,29.9 L 123.7,27.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-23 atom-20 atom-23' d='M 123.7,27.6 L 132.4,25.2' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 36.9,115.6 L 37.7,114.8 L 38.1,115.0' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 70.2,86.7 L 70.0,85.7 L 70.8,84.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 85.3,71.8 L 86.2,71.1 L 87.2,71.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 105.8,77.4 L 106.9,77.8 L 107.1,78.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 111.2,98.0 L 111.4,99.1 L 110.6,99.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 99.6,133.9 L 99.8,134.9 L 99.5,135.2' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 35.4,31.7 L 36.5,31.4 L 36.8,31.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 78.2,21.2 L 78.5,20.1 L 79.6,19.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 98.5,14.8 L 99.5,14.5 L 100.3,15.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 109.6,49.9 L 109.3,50.9 L 108.2,51.2' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 89.3,56.3 L 88.3,56.6 L 87.5,55.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-0' d='M 30.3 93.6
Q 30.3 92.1, 31.0 91.2
Q 31.8 90.4, 33.1 90.4
Q 34.5 90.4, 35.2 91.2
Q 36.0 92.1, 36.0 93.6
Q 36.0 95.1, 35.2 95.9
Q 34.5 96.8, 33.1 96.8
Q 31.8 96.8, 31.0 95.9
Q 30.3 95.1, 30.3 93.6
M 33.1 96.1
Q 34.1 96.1, 34.6 95.4
Q 35.1 94.8, 35.1 93.6
Q 35.1 92.3, 34.6 91.7
Q 34.1 91.1, 33.1 91.1
Q 32.2 91.1, 31.7 91.7
Q 31.2 92.3, 31.2 93.6
Q 31.2 94.8, 31.7 95.4
Q 32.2 96.1, 33.1 96.1
' fill='#FF0000'/>
<path class='atom-3' d='M 55.6 121.5
Q 55.6 120.1, 56.3 119.2
Q 57.0 118.4, 58.4 118.4
Q 59.8 118.4, 60.5 119.2
Q 61.2 120.1, 61.2 121.5
Q 61.2 123.0, 60.5 123.9
Q 59.7 124.7, 58.4 124.7
Q 57.0 124.7, 56.3 123.9
Q 55.6 123.0, 55.6 121.5
M 58.4 124.0
Q 59.3 124.0, 59.8 123.4
Q 60.3 122.8, 60.3 121.5
Q 60.3 120.3, 59.8 119.7
Q 59.3 119.1, 58.4 119.1
Q 57.4 119.1, 56.9 119.7
Q 56.4 120.3, 56.4 121.5
Q 56.4 122.8, 56.9 123.4
Q 57.4 124.0, 58.4 124.0
' fill='#FF0000'/>
<path class='atom-11' d='M 117.7 141.7
Q 117.7 140.2, 118.4 139.3
Q 119.1 138.5, 120.5 138.5
Q 121.9 138.5, 122.6 139.3
Q 123.3 140.2, 123.3 141.7
Q 123.3 143.1, 122.6 144.0
Q 121.9 144.8, 120.5 144.8
Q 119.2 144.8, 118.4 144.0
Q 117.7 143.2, 117.7 141.7
M 120.5 144.2
Q 121.5 144.2, 122.0 143.5
Q 122.5 142.9, 122.5 141.7
Q 122.5 140.4, 122.0 139.8
Q 121.5 139.2, 120.5 139.2
Q 119.6 139.2, 119.1 139.8
Q 118.6 140.4, 118.6 141.7
Q 118.6 142.9, 119.1 143.5
Q 119.6 144.2, 120.5 144.2
' fill='#FF0000'/>
<path class='atom-12' d='M 74.7 146.5
L 75.6 146.5
L 75.6 149.1
L 78.7 149.1
L 78.7 146.5
L 79.5 146.5
L 79.5 152.6
L 78.7 152.6
L 78.7 149.8
L 75.6 149.8
L 75.6 152.6
L 74.7 152.6
L 74.7 146.5
' fill='#FF0000'/>
<path class='atom-12' d='M 80.8 149.5
Q 80.8 148.1, 81.5 147.2
Q 82.3 146.4, 83.6 146.4
Q 85.0 146.4, 85.7 147.2
Q 86.5 148.1, 86.5 149.5
Q 86.5 151.0, 85.7 151.9
Q 85.0 152.7, 83.6 152.7
Q 82.3 152.7, 81.5 151.9
Q 80.8 151.0, 80.8 149.5
M 83.6 152.0
Q 84.6 152.0, 85.1 151.4
Q 85.6 150.8, 85.6 149.5
Q 85.6 148.3, 85.1 147.7
Q 84.6 147.1, 83.6 147.1
Q 82.7 147.1, 82.2 147.7
Q 81.7 148.3, 81.7 149.5
Q 81.7 150.8, 82.2 151.4
Q 82.7 152.0, 83.6 152.0
' fill='#FF0000'/>
<path class='atom-15' d='M 39.3 10.4
Q 39.3 8.9, 40.0 8.1
Q 40.7 7.3, 42.1 7.3
Q 43.5 7.3, 44.2 8.1
Q 44.9 8.9, 44.9 10.4
Q 44.9 11.9, 44.2 12.8
Q 43.4 13.6, 42.1 13.6
Q 40.7 13.6, 40.0 12.8
Q 39.3 11.9, 39.3 10.4
M 42.1 12.9
Q 43.0 12.9, 43.5 12.3
Q 44.1 11.6, 44.1 10.4
Q 44.1 9.2, 43.5 8.6
Q 43.0 8.0, 42.1 8.0
Q 41.2 8.0, 40.6 8.6
Q 40.1 9.2, 40.1 10.4
Q 40.1 11.7, 40.6 12.3
Q 41.2 12.9, 42.1 12.9
' fill='#FF0000'/>
<path class='atom-16' d='M 50.5 43.7
L 52.5 47.0
Q 52.7 47.3, 53.0 47.9
Q 53.4 48.5, 53.4 48.5
L 53.4 43.7
L 54.2 43.7
L 54.2 49.9
L 53.3 49.9
L 51.2 46.3
Q 50.9 45.9, 50.7 45.4
Q 50.4 44.9, 50.3 44.8
L 50.3 49.9
L 49.5 49.9
L 49.5 43.7
L 50.5 43.7
' fill='#0000FF'/>
<path class='atom-16' d='M 49.4 50.5
L 50.3 50.5
L 50.3 53.1
L 53.4 53.1
L 53.4 50.5
L 54.3 50.5
L 54.3 56.7
L 53.4 56.7
L 53.4 53.8
L 50.3 53.8
L 50.3 56.7
L 49.4 56.7
L 49.4 50.5
' fill='#0000FF'/>
<path class='atom-23' d='M 133.1 24.3
Q 133.1 22.8, 133.9 22.0
Q 134.6 21.2, 136.0 21.2
Q 137.3 21.2, 138.0 22.0
Q 138.8 22.8, 138.8 24.3
Q 138.8 25.8, 138.0 26.6
Q 137.3 27.5, 136.0 27.5
Q 134.6 27.5, 133.9 26.6
Q 133.1 25.8, 133.1 24.3
M 136.0 26.8
Q 136.9 26.8, 137.4 26.2
Q 137.9 25.5, 137.9 24.3
Q 137.9 23.1, 137.4 22.5
Q 136.9 21.9, 136.0 21.9
Q 135.0 21.9, 134.5 22.5
Q 134.0 23.1, 134.0 24.3
Q 134.0 25.5, 134.5 26.2
Q 135.0 26.8, 136.0 26.8
' fill='#FF0000'/>
<path class='atom-23' d='M 139.7 21.2
L 140.6 21.2
L 140.6 23.8
L 143.7 23.8
L 143.7 21.2
L 144.6 21.2
L 144.6 27.4
L 143.7 27.4
L 143.7 24.5
L 140.6 24.5
L 140.6 27.4
L 139.7 27.4
L 139.7 21.2
' fill='#FF0000'/>
</svg>
</span>
</div>


## Developing
This repo uses pre-commit, so after cloning run `pip install -r requirements.txt` and `pre-commit install` prior to committing.
