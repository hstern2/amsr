# AMSR
**A**nother **M**olecular **S**tring **R**epresentation,
inspired by

- [H. Hiz, "A Linearization of Chemical Graphs," *J. Chem. Doc.* **4**, 173-180 (1964)](https://pubs.acs.org/doi/10.1021/c160014a015)
- [SMILES](https://pubs.acs.org/doi/10.1021/ci00057a005)
- [PATTY](https://pubs.acs.org/doi/10.1021/ci00015a015)
- [DeepSMILES](https://github.com/baoilleach/deepsmiles)
- [SELFIES](https://github.com/aspuru-guzik-group/selfies)

## Demo
[AMSR demo](http://155.138.219.61:8000)

## Installing
```
pip install amsr
```

For GPU:

```
conda install pytorch torchvision torchaudio cudatoolkit=11.8 -c pytorch -c nvidia
```

## Usage
```py
import amsr

amsr.ToMol("CNcncc5cNcN6C.oC.o") # caffeine
```
![caffeine](https://user-images.githubusercontent.com/19351218/151638119-b1439d47-5e5a-417e-9254-c34568e2f3d1.png)

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
<path class='bond-0 atom-0 atom-1' d='M 63.3,104.4 L 105.6,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 105.6,80.0 L 123.9,69.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 123.9,69.4 L 142.2,58.9' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-1 atom-3' d='M 105.6,80.0 L 115.3,96.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-1 atom-3' d='M 115.3,96.9 L 125.1,113.8' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-1 atom-4' d='M 105.6,80.0 L 95.8,63.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-1 atom-4' d='M 95.8,63.1 L 86.1,46.2' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-0 atom-5' d='M 63.3,104.4 L 43.4,92.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-0 atom-5' d='M 43.4,92.9 L 23.5,81.4' style='fill:none;fill-rule:evenodd;stroke:#00CC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 65.4,103.2 L 63.3,104.4 L 62.3,103.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-2' d='M 143.8 48.7
L 152.0 48.7
L 152.0 50.2
L 145.6 50.2
L 145.6 54.4
L 151.3 54.4
L 151.3 56.0
L 145.6 56.0
L 145.6 62.5
L 143.8 62.5
L 143.8 48.7
' fill='#33CCCC'/>
<path class='atom-3' d='M 125.9 115.4
L 134.1 115.4
L 134.1 117.0
L 127.7 117.0
L 127.7 121.2
L 133.4 121.2
L 133.4 122.8
L 127.7 122.8
L 127.7 129.2
L 125.9 129.2
L 125.9 115.4
' fill='#33CCCC'/>
<path class='atom-4' d='M 77.0 30.8
L 85.3 30.8
L 85.3 32.4
L 78.9 32.4
L 78.9 36.6
L 84.6 36.6
L 84.6 38.2
L 78.9 38.2
L 78.9 44.6
L 77.0 44.6
L 77.0 30.8
' fill='#33CCCC'/>
<path class='atom-5' d='M 8.0 80.5
Q 8.0 77.0, 9.6 75.2
Q 11.2 73.4, 14.3 73.4
Q 17.1 73.4, 18.7 75.4
L 17.4 76.5
Q 16.3 75.0, 14.3 75.0
Q 12.2 75.0, 11.1 76.4
Q 10.0 77.8, 10.0 80.5
Q 10.0 83.2, 11.1 84.6
Q 12.3 86.0, 14.5 86.0
Q 16.0 86.0, 17.8 85.1
L 18.4 86.6
Q 17.6 87.0, 16.5 87.3
Q 15.4 87.6, 14.2 87.6
Q 11.2 87.6, 9.6 85.8
Q 8.0 83.9, 8.0 80.5
' fill='#00CC00'/>
<path class='atom-5' d='M 20.1 72.6
L 21.9 72.6
L 21.9 87.4
L 20.1 87.4
L 20.1 72.6
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
<path class='bond-0 atom-0 atom-1' d='M 8.0,142.4 L 80.0,100.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 80.0,100.8 L 152.0,142.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-1 atom-3' d='M 80.0,100.8 L 80.0,17.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
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
<path class='bond-0 atom-0 atom-1' d='M 56.0,52.3 L 8.0,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 56.0,52.3 L 104.0,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 104.0,80.0 L 131.7,32.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-2 atom-4' d='M 104.0,80.0 L 76.3,128.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-2 atom-5' d='M 104.0,80.0 L 152.0,107.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 53.6,53.7 L 56.0,52.3 L 58.4,53.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
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
<path class='bond-0 atom-0 atom-1' d='M 140.8,72.6 L 32.2,9.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 32.2,9.9 L 32.2,63.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 32.2,63.4 L 32.2,116.9' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-0' d='M 49.3,125.5 L 95.1,99.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-0' d='M 95.1,99.0 L 140.8,72.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 135.4,69.5 L 140.8,72.6 L 138.5,73.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 37.6,13.0 L 32.2,9.9 L 32.2,12.6' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-2' d='M 19.2 135.4
Q 19.2 128.6, 22.5 124.8
Q 25.9 121.0, 32.2 121.0
Q 38.4 121.0, 41.8 124.8
Q 45.2 128.6, 45.2 135.4
Q 45.2 142.3, 41.8 146.2
Q 38.4 150.1, 32.2 150.1
Q 25.9 150.1, 22.5 146.2
Q 19.2 142.4, 19.2 135.4
M 32.2 146.9
Q 36.5 146.9, 38.8 144.0
Q 41.2 141.1, 41.2 135.4
Q 41.2 129.9, 38.8 127.1
Q 36.5 124.2, 32.2 124.2
Q 27.8 124.2, 25.5 127.0
Q 23.2 129.8, 23.2 135.4
Q 23.2 141.2, 25.5 144.0
Q 27.8 146.9, 32.2 146.9
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
<path class='bond-0 atom-0 atom-1' d='M 152.0,80.0 L 116.0,142.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 116.0,142.4 L 44.0,142.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 44.0,142.4 L 8.0,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 8.0,80.0 L 44.0,17.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 44.0,17.6 L 116.0,17.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-0' d='M 116.0,17.6 L 152.0,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 150.2,83.1 L 152.0,80.0 L 150.2,76.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 117.8,139.2 L 116.0,142.4 L 112.4,142.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 47.6,142.4 L 44.0,142.4 L 42.2,139.2' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 9.8,83.1 L 8.0,80.0 L 9.8,76.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 42.2,20.8 L 44.0,17.6 L 47.6,17.6' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 112.4,17.6 L 116.0,17.6 L 117.8,20.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
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
<path class='bond-0 atom-0 atom-1' d='M 152.0,80.0 L 142.4,116.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 142.4,116.0 L 116.0,142.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 116.0,142.4 L 80.0,152.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 80.0,152.0 L 44.0,142.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 44.0,142.4 L 17.6,116.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-6' d='M 17.6,116.0 L 8.0,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-6 atom-7' d='M 8.0,80.0 L 17.6,44.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-7 atom-8' d='M 17.6,44.0 L 44.0,17.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-8 atom-9' d='M 44.0,17.6 L 80.0,8.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-9 atom-10' d='M 80.0,8.0 L 116.0,17.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-10 atom-11' d='M 116.0,17.6 L 142.4,44.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-11 atom-0' d='M 142.4,44.0 L 152.0,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 151.5,81.8 L 152.0,80.0 L 151.5,78.2' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 142.8,114.2 L 142.4,116.0 L 141.0,117.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 117.3,141.0 L 116.0,142.4 L 114.2,142.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 81.8,151.5 L 80.0,152.0 L 78.2,151.5' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 45.8,142.8 L 44.0,142.4 L 42.7,141.0' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 19.0,117.3 L 17.6,116.0 L 17.2,114.2' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 8.5,81.8 L 8.0,80.0 L 8.5,78.2' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 17.2,45.8 L 17.6,44.0 L 19.0,42.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 42.7,19.0 L 44.0,17.6 L 45.8,17.2' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 78.2,8.5 L 80.0,8.0 L 81.8,8.5' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 114.2,17.2 L 116.0,17.6 L 117.3,19.0' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 141.0,42.7 L 142.4,44.0 L 142.8,45.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
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
<path class='bond-0 atom-0 atom-1' d='M 9.8,70.2 L 64.9,70.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 64.9,70.2 L 120.0,70.2' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 9.8,89.3 L 64.9,89.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 64.9,89.3 L 120.0,89.3' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='atom-1' d='M 124.2 79.9
Q 124.2 73.1, 127.6 69.3
Q 130.9 65.5, 137.2 65.5
Q 143.5 65.5, 146.9 69.3
Q 150.2 73.1, 150.2 79.9
Q 150.2 86.7, 146.8 90.7
Q 143.4 94.5, 137.2 94.5
Q 131.0 94.5, 127.6 90.7
Q 124.2 86.8, 124.2 79.9
M 137.2 91.3
Q 141.5 91.3, 143.9 88.5
Q 146.2 85.5, 146.2 79.9
Q 146.2 74.3, 143.9 71.5
Q 141.5 68.7, 137.2 68.7
Q 132.9 68.7, 130.5 71.5
Q 128.2 74.3, 128.2 79.9
Q 128.2 85.6, 130.5 88.5
Q 132.9 91.3, 137.2 91.3
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
<path class='bond-0 atom-0 atom-1' d='M 152.0,80.0 L 116.0,142.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 139.5,80.0 L 109.8,131.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 116.0,142.4 L 44.0,142.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 44.0,142.4 L 8.0,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 50.2,131.6 L 20.5,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 8.0,80.0 L 44.0,17.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 44.0,17.6 L 116.0,17.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 50.2,28.4 L 109.8,28.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-0' d='M 116.0,17.6 L 152.0,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 150.2,83.1 L 152.0,80.0 L 150.2,76.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 117.8,139.2 L 116.0,142.4 L 112.4,142.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 47.6,142.4 L 44.0,142.4 L 42.2,139.2' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 9.8,83.1 L 8.0,80.0 L 9.8,76.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 42.2,20.8 L 44.0,17.6 L 47.6,17.6' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 112.4,17.6 L 116.0,17.6 L 117.8,20.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
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
<path class='bond-0 atom-0 atom-1' d='M 8.0,91.3 L 76.3,51.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 76.3,51.9 L 104.1,67.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 104.1,67.9 L 131.8,83.9' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 76.3,65.6 L 98.2,78.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 98.2,78.2 L 125.9,94.2' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 72.9,53.9 L 76.3,51.9 L 77.7,52.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-2' d='M 131.5 96.5
Q 131.5 91.2, 134.1 88.2
Q 136.8 85.2, 141.7 85.2
Q 146.7 85.2, 149.3 88.2
Q 152.0 91.2, 152.0 96.5
Q 152.0 102.0, 149.3 105.1
Q 146.6 108.1, 141.7 108.1
Q 136.8 108.1, 134.1 105.1
Q 131.5 102.0, 131.5 96.5
M 141.7 105.6
Q 145.1 105.6, 147.0 103.3
Q 148.8 101.0, 148.8 96.5
Q 148.8 92.1, 147.0 89.9
Q 145.1 87.7, 141.7 87.7
Q 138.3 87.7, 136.5 89.9
Q 134.6 92.1, 134.6 96.5
Q 134.6 101.0, 136.5 103.3
Q 138.3 105.6, 141.7 105.6
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
<path class='bond-0 atom-0 atom-1' d='M 142.8,74.0 L 94.8,8.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 127.6,74.0 L 90.1,22.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 94.8,8.0 L 17.2,33.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 17.2,33.2 L 17.2,114.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 29.5,42.1 L 29.5,105.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 17.2,114.8 L 49.4,125.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 49.4,125.2 L 81.5,135.6' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-0' d='M 105.3,125.6 L 124.0,99.8' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-0' d='M 124.0,99.8 L 142.8,74.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 140.4,70.7 L 142.8,74.0 L 141.8,75.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 97.2,11.3 L 94.8,8.0 L 90.9,9.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 21.1,31.9 L 17.2,33.2 L 17.2,37.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 17.2,110.7 L 17.2,114.8 L 18.9,115.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-4' d='M 84.2 140.0
Q 84.2 134.5, 87.0 131.4
Q 89.7 128.3, 94.8 128.3
Q 99.9 128.3, 102.7 131.4
Q 105.4 134.5, 105.4 140.0
Q 105.4 145.6, 102.6 148.8
Q 99.9 152.0, 94.8 152.0
Q 89.7 152.0, 87.0 148.8
Q 84.2 145.7, 84.2 140.0
M 94.8 149.4
Q 98.3 149.4, 100.2 147.0
Q 102.2 144.7, 102.2 140.0
Q 102.2 135.5, 100.2 133.2
Q 98.3 130.9, 94.8 130.9
Q 91.3 130.9, 89.4 133.2
Q 87.5 135.5, 87.5 140.0
Q 87.5 144.7, 89.4 147.0
Q 91.3 149.4, 94.8 149.4
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
<path class='bond-0 atom-0 atom-1' d='M 133.5,64.2 L 92.6,8.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 120.6,64.2 L 88.6,20.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 92.6,8.0 L 26.5,29.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 26.5,29.5 L 26.5,99.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 36.9,37.1 L 36.9,91.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 26.5,99.0 L 54.7,108.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 54.7,108.2 L 82.9,117.3' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-0' d='M 101.4,108.4 L 117.5,86.3' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-0' d='M 117.5,86.3 L 133.5,64.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 131.5,61.4 L 133.5,64.2 L 132.7,65.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 94.7,10.8 L 92.6,8.0 L 89.3,9.1' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 29.8,28.4 L 26.5,29.5 L 26.5,33.0' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 26.5,95.5 L 26.5,99.0 L 27.9,99.5' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-4' d='M 88.3 110.7
L 94.7 121.1
Q 95.4 122.1, 96.4 124.0
Q 97.4 125.8, 97.5 125.9
L 97.5 110.7
L 100.1 110.7
L 100.1 130.3
L 97.4 130.3
L 90.5 118.9
Q 89.7 117.6, 88.8 116.1
Q 88.0 114.5, 87.7 114.1
L 87.7 130.3
L 85.2 130.3
L 85.2 110.7
L 88.3 110.7
' fill='#0000FF'/>
<path class='atom-4' d='M 84.9 132.3
L 87.6 132.3
L 87.6 140.7
L 97.7 140.7
L 97.7 132.3
L 100.3 132.3
L 100.3 152.0
L 97.7 152.0
L 97.7 142.9
L 87.6 142.9
L 87.6 152.0
L 84.9 152.0
L 84.9 132.3
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
<path class='bond-0 atom-0 atom-1' d='M 33.9,115.6 L 8.0,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 36.4,107.9 L 16.2,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 8.0,80.0 L 18.1,66.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 18.1,66.2 L 28.1,52.3' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 41.1,46.7 L 58.4,52.4' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 58.4,52.4 L 75.7,58.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 75.7,58.0 L 75.7,102.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 69.1,62.8 L 69.1,97.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-0' d='M 75.7,102.0 L 33.9,115.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-4 atom-5' d='M 75.7,102.0 L 113.9,124.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-5 atom-6' d='M 113.9,124.0 L 152.0,102.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-5 atom-6' d='M 113.9,116.4 L 145.4,98.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-6 atom-7' d='M 152.0,102.0 L 152.0,58.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-7 atom-8' d='M 152.0,58.0 L 113.9,36.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-7 atom-8' d='M 145.4,61.8 L 113.9,43.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-8 atom-3' d='M 113.9,36.0 L 75.7,58.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 32.6,113.8 L 33.9,115.6 L 36.0,114.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 9.3,81.8 L 8.0,80.0 L 8.5,79.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 112.0,122.9 L 113.9,124.0 L 115.8,122.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 150.1,103.1 L 152.0,102.0 L 152.0,99.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 152.0,60.2 L 152.0,58.0 L 150.1,56.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 115.8,37.1 L 113.9,36.0 L 112.0,37.1' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-2' d='M 28.2 44.4
Q 28.2 41.4, 29.6 39.8
Q 31.1 38.1, 33.9 38.1
Q 36.6 38.1, 38.1 39.8
Q 39.6 41.4, 39.6 44.4
Q 39.6 47.4, 38.1 49.2
Q 36.6 50.9, 33.9 50.9
Q 31.1 50.9, 29.6 49.2
Q 28.2 47.5, 28.2 44.4
M 33.9 49.5
Q 35.8 49.5, 36.8 48.2
Q 37.8 46.9, 37.8 44.4
Q 37.8 42.0, 36.8 40.7
Q 35.8 39.5, 33.9 39.5
Q 32.0 39.5, 30.9 40.7
Q 29.9 42.0, 29.9 44.4
Q 29.9 46.9, 30.9 48.2
Q 32.0 49.5, 33.9 49.5
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
<path class='bond-0 atom-0 atom-1' d='M 78.7,58.8 L 38.4,45.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 72.3,63.4 L 40.8,53.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 38.4,45.7 L 28.7,59.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 28.7,59.1 L 18.9,72.5' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 19.1,87.6 L 28.7,101.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 28.7,101.0 L 38.4,114.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 38.4,114.3 L 78.7,101.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 40.8,106.8 L 72.3,96.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-0' d='M 78.7,101.2 L 78.7,58.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-4 atom-5' d='M 78.7,101.2 L 115.3,122.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-5 atom-6' d='M 115.3,122.3 L 152.0,101.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-5 atom-6' d='M 115.3,115.0 L 145.6,97.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-6 atom-7' d='M 152.0,101.2 L 152.0,58.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-7 atom-8' d='M 152.0,58.8 L 115.3,37.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-7 atom-8' d='M 145.6,62.5 L 115.3,45.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-8 atom-0' d='M 115.3,37.7 L 78.7,58.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 40.4,46.4 L 38.4,45.7 L 37.9,46.4' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 37.9,113.6 L 38.4,114.3 L 40.4,113.6' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 113.5,121.3 L 115.3,122.3 L 117.2,121.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 150.2,102.2 L 152.0,101.2 L 152.0,99.1' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 152.0,60.9 L 152.0,58.8 L 150.2,57.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 117.2,38.7 L 115.3,37.7 L 113.5,38.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-2' d='M 8.0 80.0
Q 8.0 77.2, 9.4 75.5
Q 10.8 73.9, 13.5 73.9
Q 16.2 73.9, 17.6 75.5
Q 19.0 77.2, 19.0 80.0
Q 19.0 82.9, 17.6 84.6
Q 16.1 86.2, 13.5 86.2
Q 10.9 86.2, 9.4 84.6
Q 8.0 83.0, 8.0 80.0
M 13.5 84.9
Q 15.3 84.9, 16.3 83.7
Q 17.3 82.4, 17.3 80.0
Q 17.3 77.7, 16.3 76.5
Q 15.3 75.3, 13.5 75.3
Q 11.7 75.3, 10.7 76.5
Q 9.7 77.7, 9.7 80.0
Q 9.7 82.5, 10.7 83.7
Q 11.7 84.9, 13.5 84.9
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
<path class='bond-0 atom-0 atom-1' d='M 150.4,80.0 L 92.9,80.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 92.9,80.0 L 35.4,80.0' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 150.4,99.5 L 92.9,99.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 92.9,99.5 L 35.4,99.5' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 150.4,60.5 L 92.9,60.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 92.9,60.5 L 35.4,60.5' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='atom-1' d='M 14.1 65.8
L 23.4 80.8
Q 24.3 82.3, 25.8 85.0
Q 27.3 87.7, 27.3 87.8
L 27.3 65.8
L 31.1 65.8
L 31.1 94.2
L 27.2 94.2
L 17.3 77.8
Q 16.1 75.8, 14.9 73.6
Q 13.7 71.4, 13.3 70.8
L 13.3 94.2
L 9.6 94.2
L 9.6 65.8
L 14.1 65.8
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
<path class='bond-0 atom-0 atom-1' d='M 133.3,84.6 L 111.2,84.6' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 111.2,84.6 L 89.1,84.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 133.3,75.1 L 111.2,75.1' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 111.2,75.1 L 89.1,75.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 70.9,84.6 L 48.8,84.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 48.8,84.6 L 26.7,84.6' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 70.9,75.1 L 48.8,75.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 48.8,75.1 L 26.7,75.1' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='atom-0' d='M 135.4 79.9
Q 135.4 75.6, 137.6 73.2
Q 139.7 70.7, 143.7 70.7
Q 147.7 70.7, 149.9 73.2
Q 152.0 75.6, 152.0 79.9
Q 152.0 84.3, 149.8 86.8
Q 147.7 89.3, 143.7 89.3
Q 139.7 89.3, 137.6 86.8
Q 135.4 84.3, 135.4 79.9
M 143.7 87.2
Q 146.5 87.2, 147.9 85.4
Q 149.5 83.5, 149.5 79.9
Q 149.5 76.4, 147.9 74.6
Q 146.5 72.8, 143.7 72.8
Q 141.0 72.8, 139.5 74.6
Q 138.0 76.3, 138.0 79.9
Q 138.0 83.6, 139.5 85.4
Q 141.0 87.2, 143.7 87.2
' fill='#FF0000'/>
<path class='atom-1' d='M 73.0 79.9
Q 73.0 75.5, 75.1 73.1
Q 77.2 70.7, 81.2 70.7
Q 85.0 70.7, 87.0 73.4
L 85.3 74.7
Q 83.8 72.8, 81.2 72.8
Q 78.5 72.8, 77.1 74.7
Q 75.6 76.5, 75.6 79.9
Q 75.6 83.5, 77.1 85.3
Q 78.6 87.2, 81.5 87.2
Q 83.5 87.2, 85.8 86.0
L 86.6 87.9
Q 85.6 88.5, 84.2 88.9
Q 82.8 89.2, 81.2 89.2
Q 77.2 89.2, 75.1 86.8
Q 73.0 84.4, 73.0 79.9
' fill='#000000'/>
<path class='atom-2' d='M 8.0 79.9
Q 8.0 75.6, 10.1 73.2
Q 12.3 70.7, 16.3 70.7
Q 20.3 70.7, 22.4 73.2
Q 24.6 75.6, 24.6 79.9
Q 24.6 84.3, 22.4 86.8
Q 20.2 89.3, 16.3 89.3
Q 12.3 89.3, 10.1 86.8
Q 8.0 84.3, 8.0 79.9
M 16.3 87.2
Q 19.0 87.2, 20.5 85.4
Q 22.0 83.5, 22.0 79.9
Q 22.0 76.4, 20.5 74.6
Q 19.0 72.8, 16.3 72.8
Q 13.5 72.8, 12.0 74.6
Q 10.5 76.3, 10.5 79.9
Q 10.5 83.6, 12.0 85.4
Q 13.5 87.2, 16.3 87.2
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
<path class='bond-0 atom-0 atom-1' d='M 8.0,106.7 L 38.9,88.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 38.9,88.8 L 69.7,71.0' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 90.3,71.0 L 121.1,88.8' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 121.1,88.8 L 152.0,106.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='atom-1' d='M 73.3 73.2
Q 73.6 73.3, 74.7 73.7
Q 75.8 74.2, 77.0 74.5
Q 78.2 74.8, 79.4 74.8
Q 81.7 74.8, 83.0 73.7
Q 84.3 72.6, 84.3 70.7
Q 84.3 69.4, 83.6 68.6
Q 83.0 67.8, 82.0 67.4
Q 81.0 66.9, 79.3 66.4
Q 77.2 65.8, 75.9 65.2
Q 74.7 64.6, 73.8 63.4
Q 72.9 62.1, 72.9 60.0
Q 72.9 57.0, 74.9 55.2
Q 77.0 53.3, 81.0 53.3
Q 83.7 53.3, 86.8 54.6
L 86.0 57.2
Q 83.2 56.0, 81.1 56.0
Q 78.8 56.0, 77.5 57.0
Q 76.2 57.9, 76.3 59.6
Q 76.3 60.8, 76.9 61.6
Q 77.6 62.4, 78.5 62.8
Q 79.5 63.2, 81.1 63.7
Q 83.2 64.4, 84.5 65.1
Q 85.7 65.7, 86.6 67.1
Q 87.5 68.4, 87.5 70.7
Q 87.5 74.0, 85.4 75.7
Q 83.2 77.5, 79.6 77.5
Q 77.5 77.5, 75.9 77.0
Q 74.3 76.6, 72.5 75.8
L 73.3 73.2
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
<path class='bond-0 atom-0 atom-1' d='M 8.0,148.3 L 38.9,130.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 38.9,130.5 L 69.7,112.7' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 86.2,92.3 L 86.2,65.4' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 86.2,65.4 L 86.2,38.6' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 73.8,92.3 L 73.8,65.4' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 73.8,65.4 L 73.8,38.6' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-1 atom-3' d='M 90.3,112.7 L 121.1,130.5' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-1 atom-3' d='M 121.1,130.5 L 152.0,148.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='atom-1' d='M 73.3 114.8
Q 73.6 114.9, 74.7 115.4
Q 75.8 115.8, 77.0 116.1
Q 78.2 116.4, 79.4 116.4
Q 81.7 116.4, 83.0 115.4
Q 84.3 114.3, 84.3 112.4
Q 84.3 111.1, 83.6 110.3
Q 83.0 109.5, 82.0 109.0
Q 81.0 108.6, 79.3 108.1
Q 77.2 107.5, 75.9 106.9
Q 74.7 106.3, 73.8 105.0
Q 72.9 103.7, 72.9 101.6
Q 72.9 98.7, 74.9 96.8
Q 77.0 95.0, 81.0 95.0
Q 83.7 95.0, 86.8 96.3
L 86.0 98.9
Q 83.2 97.7, 81.1 97.7
Q 78.8 97.7, 77.5 98.7
Q 76.2 99.6, 76.3 101.2
Q 76.3 102.5, 76.9 103.2
Q 77.6 104.0, 78.5 104.4
Q 79.5 104.9, 81.1 105.4
Q 83.2 106.0, 84.5 106.7
Q 85.7 107.4, 86.6 108.7
Q 87.5 110.1, 87.5 112.4
Q 87.5 115.6, 85.4 117.4
Q 83.2 119.1, 79.6 119.1
Q 77.5 119.1, 75.9 118.6
Q 74.3 118.2, 72.5 117.4
L 73.3 114.8
' fill='#CCCC00'/>
<path class='atom-2' d='M 69.2 23.7
Q 69.2 18.0, 72.0 14.9
Q 74.8 11.7, 80.0 11.7
Q 85.2 11.7, 88.0 14.9
Q 90.8 18.0, 90.8 23.7
Q 90.8 29.4, 88.0 32.6
Q 85.2 35.9, 80.0 35.9
Q 74.8 35.9, 72.0 32.6
Q 69.2 29.4, 69.2 23.7
M 80.0 33.2
Q 83.6 33.2, 85.5 30.8
Q 87.5 28.4, 87.5 23.7
Q 87.5 19.0, 85.5 16.7
Q 83.6 14.4, 80.0 14.4
Q 76.4 14.4, 74.4 16.7
Q 72.5 19.0, 72.5 23.7
Q 72.5 28.4, 74.4 30.8
Q 76.4 33.2, 80.0 33.2
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
<path class='bond-0 atom-0 atom-1' d='M 72.2,84.5 L 52.5,95.9' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 52.5,95.9 L 32.8,107.3' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 80.0,91.5 L 80.0,111.8' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 80.0,111.8 L 80.0,132.1' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 87.8,84.5 L 107.5,95.9' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 107.5,95.9 L 127.2,107.3' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-0 atom-4' d='M 80.0,69.0 L 80.0,48.5' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-0 atom-4' d='M 80.0,48.5 L 80.0,27.9' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-0 atom-5' d='M 87.8,75.5 L 107.5,64.1' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-0 atom-5' d='M 107.5,64.1 L 127.2,52.7' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-0 atom-6' d='M 72.2,75.5 L 52.5,64.1' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-0 atom-6' d='M 52.5,64.1 L 32.8,52.7' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='atom-0' d='M 75.0 86.1
Q 75.2 86.2, 76.0 86.6
Q 76.8 86.9, 77.7 87.1
Q 78.7 87.3, 79.6 87.3
Q 81.3 87.3, 82.2 86.5
Q 83.2 85.7, 83.2 84.3
Q 83.2 83.3, 82.7 82.7
Q 82.2 82.1, 81.5 81.7
Q 80.7 81.4, 79.5 81.0
Q 77.9 80.6, 76.9 80.1
Q 76.0 79.6, 75.3 78.7
Q 74.7 77.7, 74.7 76.1
Q 74.7 73.9, 76.2 72.5
Q 77.7 71.1, 80.7 71.1
Q 82.8 71.1, 85.1 72.1
L 84.6 74.0
Q 82.4 73.1, 80.8 73.1
Q 79.1 73.1, 78.1 73.9
Q 77.1 74.6, 77.2 75.8
Q 77.2 76.8, 77.7 77.4
Q 78.2 77.9, 78.9 78.3
Q 79.6 78.6, 80.8 79.0
Q 82.4 79.5, 83.4 80.0
Q 84.3 80.5, 85.0 81.5
Q 85.7 82.5, 85.7 84.3
Q 85.7 86.7, 84.1 88.1
Q 82.4 89.4, 79.7 89.4
Q 78.1 89.4, 76.9 89.0
Q 75.7 88.7, 74.3 88.1
L 75.0 86.1
' fill='#CCCC00'/>
<path class='atom-1' d='M 20.1 102.6
L 30.7 102.6
L 30.7 104.6
L 22.5 104.6
L 22.5 110.1
L 29.8 110.1
L 29.8 112.1
L 22.5 112.1
L 22.5 120.5
L 20.1 120.5
L 20.1 102.6
' fill='#33CCCC'/>
<path class='atom-2' d='M 74.7 134.1
L 85.3 134.1
L 85.3 136.2
L 77.1 136.2
L 77.1 141.6
L 84.4 141.6
L 84.4 143.7
L 77.1 143.7
L 77.1 152.0
L 74.7 152.0
L 74.7 134.1
' fill='#33CCCC'/>
<path class='atom-3' d='M 129.3 102.6
L 139.9 102.6
L 139.9 104.6
L 131.7 104.6
L 131.7 110.1
L 139.0 110.1
L 139.0 112.1
L 131.7 112.1
L 131.7 120.5
L 129.3 120.5
L 129.3 102.6
' fill='#33CCCC'/>
<path class='atom-4' d='M 74.7 8.0
L 85.3 8.0
L 85.3 10.0
L 77.1 10.0
L 77.1 15.5
L 84.4 15.5
L 84.4 17.5
L 77.1 17.5
L 77.1 25.9
L 74.7 25.9
L 74.7 8.0
' fill='#33CCCC'/>
<path class='atom-5' d='M 129.3 39.5
L 139.9 39.5
L 139.9 41.6
L 131.7 41.6
L 131.7 47.0
L 139.0 47.0
L 139.0 49.1
L 131.7 49.1
L 131.7 57.4
L 129.3 57.4
L 129.3 39.5
' fill='#33CCCC'/>
<path class='atom-6' d='M 20.1 39.5
L 30.7 39.5
L 30.7 41.6
L 22.5 41.6
L 22.5 47.0
L 29.8 47.0
L 29.8 49.1
L 22.5 49.1
L 22.5 57.4
L 20.1 57.4
L 20.1 39.5
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
[Mg++];S!!\:ooO-O- | [magnesium sulfate](https://en.wikipedia.org/wiki/Magnesium_sulfate)
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
<path class='bond-0 atom-1 atom-2' d='M 44.8,80.2 L 32.0,87.6' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-1 atom-2' d='M 32.0,87.6 L 19.3,94.9' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-1 atom-2' d='M 48.1,85.9 L 35.3,93.3' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-1 atom-2' d='M 35.3,93.3 L 22.6,100.6' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-3' d='M 59.0,79.6 L 71.8,72.2' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-3' d='M 71.8,72.2 L 84.5,64.9' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-3' d='M 55.7,73.9 L 68.4,66.5' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-3' d='M 68.4,66.5 L 81.2,59.2' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-1 atom-4' d='M 56.5,87.9 L 63.0,99.1' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-1 atom-4' d='M 63.0,99.1 L 69.4,110.3' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-1 atom-5' d='M 47.5,72.2 L 41.0,61.0' style='fill:none;fill-rule:evenodd;stroke:#CCCC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-1 atom-5' d='M 41.0,61.0 L 34.4,49.7' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='atom-0' d='M 126.5 86.1
L 124.9 86.1
L 123.6 75.9
L 120.3 86.1
L 118.6 86.1
L 115.2 75.9
L 113.9 86.1
L 112.3 86.1
L 113.9 73.7
L 116.0 73.7
L 119.4 83.7
L 122.9 73.7
L 125.0 73.7
L 126.5 86.1
' fill='#000000'/>
<path class='atom-0' d='M 136.9 77.1
L 136.9 86.7
Q 136.9 90.6, 132.8 90.6
Q 130.8 90.6, 129.0 89.6
L 129.6 88.4
Q 130.5 88.9, 131.2 89.1
Q 131.9 89.3, 132.8 89.3
Q 134.1 89.3, 134.7 88.7
Q 135.2 88.1, 135.2 86.7
L 135.2 85.1
Q 134.3 86.3, 132.7 86.3
Q 130.8 86.3, 129.7 85.1
Q 128.7 83.9, 128.7 81.8
Q 128.7 79.5, 129.9 78.2
Q 131.2 76.9, 133.3 76.9
Q 134.5 76.9, 135.5 77.3
L 135.6 77.1
L 136.9 77.1
M 132.8 85.0
Q 133.7 85.0, 134.3 84.5
Q 135.0 84.0, 135.2 83.1
L 135.2 78.5
Q 134.3 78.2, 133.3 78.2
Q 131.9 78.2, 131.2 79.2
Q 130.4 80.1, 130.4 81.8
Q 130.4 83.3, 131.0 84.2
Q 131.6 85.0, 132.8 85.0
' fill='#000000'/>
<path class='atom-0' d='M 139.4 72.6
Q 139.7 71.9, 140.4 71.4
Q 141.1 71.0, 142.1 71.0
Q 143.4 71.0, 144.1 71.7
Q 144.8 72.3, 144.8 73.5
Q 144.8 74.8, 143.8 75.9
Q 142.9 77.0, 141.1 78.4
L 144.9 78.4
L 144.9 79.3
L 139.3 79.3
L 139.3 78.5
Q 140.9 77.4, 141.8 76.6
Q 142.7 75.8, 143.1 75.1
Q 143.6 74.4, 143.6 73.6
Q 143.6 72.8, 143.2 72.4
Q 142.8 71.9, 142.1 71.9
Q 141.4 71.9, 141.0 72.2
Q 140.6 72.5, 140.3 73.1
L 139.4 72.6
' fill='#000000'/>
<path class='atom-0' d='M 146.6 75.3
L 148.8 75.3
L 148.8 73.0
L 149.7 73.0
L 149.7 75.3
L 152.0 75.3
L 152.0 76.2
L 149.7 76.2
L 149.7 78.5
L 148.8 78.5
L 148.8 76.2
L 146.6 76.2
L 146.6 75.3
' fill='#000000'/>
<path class='atom-1' d='M 48.4 84.2
Q 48.5 84.2, 49.1 84.5
Q 49.7 84.7, 50.3 84.9
Q 51.0 85.0, 51.6 85.0
Q 52.8 85.0, 53.5 84.5
Q 54.1 83.9, 54.1 82.9
Q 54.1 82.2, 53.8 81.8
Q 53.5 81.3, 52.9 81.1
Q 52.4 80.9, 51.5 80.6
Q 50.4 80.3, 49.7 80.0
Q 49.1 79.7, 48.6 79.0
Q 48.1 78.3, 48.1 77.2
Q 48.1 75.6, 49.2 74.7
Q 50.3 73.7, 52.4 73.7
Q 53.8 73.7, 55.5 74.4
L 55.1 75.7
Q 53.6 75.1, 52.5 75.1
Q 51.2 75.1, 50.6 75.6
Q 49.9 76.1, 49.9 77.0
Q 49.9 77.6, 50.2 78.1
Q 50.6 78.5, 51.1 78.7
Q 51.6 78.9, 52.5 79.2
Q 53.6 79.5, 54.2 79.9
Q 54.9 80.2, 55.4 81.0
Q 55.9 81.7, 55.9 82.9
Q 55.9 84.6, 54.7 85.5
Q 53.6 86.5, 51.7 86.5
Q 50.5 86.5, 49.7 86.2
Q 48.9 86.0, 47.9 85.6
L 48.4 84.2
' fill='#CCCC00'/>
<path class='atom-2' d='M 8.0 102.0
Q 8.0 99.0, 9.5 97.3
Q 11.0 95.6, 13.7 95.6
Q 16.5 95.6, 18.0 97.3
Q 19.5 99.0, 19.5 102.0
Q 19.5 105.0, 18.0 106.7
Q 16.5 108.4, 13.7 108.4
Q 11.0 108.4, 9.5 106.7
Q 8.0 105.0, 8.0 102.0
M 13.7 107.0
Q 15.6 107.0, 16.7 105.8
Q 17.7 104.5, 17.7 102.0
Q 17.7 99.5, 16.7 98.3
Q 15.6 97.0, 13.7 97.0
Q 11.8 97.0, 10.8 98.3
Q 9.8 99.5, 9.8 102.0
Q 9.8 104.5, 10.8 105.8
Q 11.8 107.0, 13.7 107.0
' fill='#FF0000'/>
<path class='atom-3' d='M 84.3 57.9
Q 84.3 54.9, 85.8 53.2
Q 87.3 51.6, 90.0 51.6
Q 92.8 51.6, 94.3 53.2
Q 95.8 54.9, 95.8 57.9
Q 95.8 60.9, 94.3 62.7
Q 92.8 64.4, 90.0 64.4
Q 87.3 64.4, 85.8 62.7
Q 84.3 61.0, 84.3 57.9
M 90.0 63.0
Q 91.9 63.0, 93.0 61.7
Q 94.0 60.4, 94.0 57.9
Q 94.0 55.5, 93.0 54.2
Q 91.9 53.0, 90.0 53.0
Q 88.1 53.0, 87.1 54.2
Q 86.1 55.4, 86.1 57.9
Q 86.1 60.4, 87.1 61.7
Q 88.1 63.0, 90.0 63.0
' fill='#FF0000'/>
<path class='atom-4' d='M 68.2 118.1
Q 68.2 115.1, 69.7 113.4
Q 71.1 111.8, 73.9 111.8
Q 76.7 111.8, 78.2 113.4
Q 79.6 115.1, 79.6 118.1
Q 79.6 121.1, 78.1 122.9
Q 76.6 124.6, 73.9 124.6
Q 71.2 124.6, 69.7 122.9
Q 68.2 121.1, 68.2 118.1
M 73.9 123.2
Q 75.8 123.2, 76.8 121.9
Q 77.9 120.6, 77.9 118.1
Q 77.9 115.6, 76.8 114.4
Q 75.8 113.2, 73.9 113.2
Q 72.0 113.2, 71.0 114.4
Q 70.0 115.6, 70.0 118.1
Q 70.0 120.6, 71.0 121.9
Q 72.0 123.2, 73.9 123.2
' fill='#FF0000'/>
<path class='atom-4' d='M 81.4 113.8
L 85.7 113.8
L 85.7 114.7
L 81.4 114.7
L 81.4 113.8
' fill='#FF0000'/>
<path class='atom-5' d='M 24.1 41.8
Q 24.1 38.8, 25.6 37.1
Q 27.1 35.4, 29.9 35.4
Q 32.6 35.4, 34.1 37.1
Q 35.6 38.8, 35.6 41.8
Q 35.6 44.8, 34.1 46.5
Q 32.6 48.2, 29.9 48.2
Q 27.1 48.2, 25.6 46.5
Q 24.1 44.8, 24.1 41.8
M 29.9 46.8
Q 31.8 46.8, 32.8 45.6
Q 33.8 44.3, 33.8 41.8
Q 33.8 39.3, 32.8 38.1
Q 31.8 36.8, 29.9 36.8
Q 28.0 36.8, 26.9 38.1
Q 25.9 39.3, 25.9 41.8
Q 25.9 44.3, 26.9 45.6
Q 28.0 46.8, 29.9 46.8
' fill='#FF0000'/>
<path class='atom-5' d='M 37.4 37.5
L 41.6 37.5
L 41.6 38.4
L 37.4 38.4
L 37.4 37.5
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
<path class='bond-0 atom-0 atom-1' d='M 46.3,77.1 L 33.0,113.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 46.3,77.1 L 8.0,70.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 46.3,77.1 L 46.3,38.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 46.3,38.2 L 80.0,18.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 80.0,18.8 L 113.7,38.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-6' d='M 113.7,38.2 L 113.7,77.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-6 atom-7' d='M 113.7,77.1 L 152.0,70.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-6 atom-8' d='M 113.7,77.1 L 127.0,113.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-6 atom-9' d='M 113.7,77.1 L 99.6,85.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-6 atom-9' d='M 99.6,85.3 L 85.5,93.4' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-9 atom-0' d='M 74.5,93.4 L 60.4,85.3' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-9 atom-0' d='M 60.4,85.3 L 46.3,77.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 80.0,103.4 L 80.0,116.0' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 80.0,116.0 L 80.0,128.6' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 46.3,40.2 L 46.3,38.2 L 48.0,37.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 78.3,19.7 L 80.0,18.8 L 81.7,19.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 112.0,37.3 L 113.7,38.2 L 113.7,40.2' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-9' d='M 77.6 91.1
L 81.2 96.9
Q 81.5 97.5, 82.1 98.5
Q 82.7 99.6, 82.7 99.6
L 82.7 91.1
L 84.2 91.1
L 84.2 102.1
L 82.7 102.1
L 78.8 95.7
Q 78.3 95.0, 77.9 94.1
Q 77.4 93.2, 77.3 93.0
L 77.3 102.1
L 75.8 102.1
L 75.8 91.1
L 77.6 91.1
' fill='#0000FF'/>
<path class='atom-10' d='M 74.9 135.5
Q 74.9 132.9, 76.2 131.4
Q 77.6 129.9, 80.0 129.9
Q 82.4 129.9, 83.8 131.4
Q 85.1 132.9, 85.1 135.5
Q 85.1 138.2, 83.7 139.7
Q 82.4 141.2, 80.0 141.2
Q 77.6 141.2, 76.2 139.7
Q 74.9 138.2, 74.9 135.5
M 80.0 140.0
Q 81.7 140.0, 82.6 138.9
Q 83.5 137.7, 83.5 135.5
Q 83.5 133.3, 82.6 132.3
Q 81.7 131.2, 80.0 131.2
Q 78.3 131.2, 77.4 132.2
Q 76.5 133.3, 76.5 135.5
Q 76.5 137.7, 77.4 138.9
Q 78.3 140.0, 80.0 140.0
' fill='#FF0000'/>
<path class='atom-10' d='M 88.2,135.5 L 88.2,135.5 L 88.2,135.6 L 88.1,135.7 L 88.1,135.7 L 88.1,135.8 L 88.1,135.9 L 88.0,135.9 L 88.0,136.0 L 87.9,136.0 L 87.9,136.1 L 87.8,136.1 L 87.8,136.1 L 87.7,136.2 L 87.7,136.2 L 87.6,136.2 L 87.5,136.2 L 87.5,136.3 L 87.4,136.3 L 87.3,136.3 L 87.3,136.3 L 87.2,136.2 L 87.1,136.2 L 87.1,136.2 L 87.0,136.2 L 87.0,136.1 L 86.9,136.1 L 86.9,136.0 L 86.8,136.0 L 86.8,135.9 L 86.7,135.9 L 86.7,135.8 L 86.7,135.8 L 86.6,135.7 L 86.6,135.6 L 86.6,135.6 L 86.6,135.5 L 86.6,135.4 L 86.6,135.4 L 86.6,135.3 L 86.6,135.3 L 86.7,135.2 L 86.7,135.1 L 86.7,135.1 L 86.8,135.0 L 86.8,135.0 L 86.9,134.9 L 86.9,134.9 L 87.0,134.8 L 87.0,134.8 L 87.1,134.8 L 87.1,134.7 L 87.2,134.7 L 87.3,134.7 L 87.3,134.7 L 87.4,134.7 L 87.5,134.7 L 87.5,134.7 L 87.6,134.7 L 87.7,134.8 L 87.7,134.8 L 87.8,134.8 L 87.8,134.8 L 87.9,134.9 L 87.9,134.9 L 88.0,135.0 L 88.0,135.0 L 88.1,135.1 L 88.1,135.2 L 88.1,135.2 L 88.1,135.3 L 88.2,135.3 L 88.2,135.4 L 88.2,135.5 L 87.4,135.5 Z' style='fill:#000000;fill-rule:evenodd;fill-opacity:1;stroke:#000000;stroke-width:0.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
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
<path class='bond-0 atom-0 atom-1' d='M 134.7,70.0 L 115.6,81.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 115.6,81.0 L 96.6,92.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 74.4,92.0 L 55.4,81.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 55.4,81.0 L 36.4,70.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='atom-0' d='M 125.9 54.1
Q 126.4 52.9, 127.5 52.2
Q 128.6 51.5, 130.1 51.5
Q 132.1 51.5, 133.1 52.6
Q 134.2 53.6, 134.2 55.5
Q 134.2 57.4, 132.8 59.1
Q 131.4 60.9, 128.6 63.0
L 134.4 63.0
L 134.4 64.4
L 125.9 64.4
L 125.9 63.2
Q 128.2 61.5, 129.6 60.2
Q 131.1 59.0, 131.7 57.9
Q 132.4 56.7, 132.4 55.6
Q 132.4 54.3, 131.8 53.7
Q 131.2 53.0, 130.1 53.0
Q 129.1 53.0, 128.4 53.4
Q 127.8 53.8, 127.3 54.7
L 125.9 54.1
' fill='#000000'/>
<path class='atom-0' d='M 136.9 54.8
L 139.5 54.8
L 139.5 63.0
L 149.4 63.0
L 149.4 54.8
L 152.0 54.8
L 152.0 74.0
L 149.4 74.0
L 149.4 65.1
L 139.5 65.1
L 139.5 74.0
L 136.9 74.0
L 136.9 54.8
' fill='#000000'/>
<path class='atom-1' d='M 76.7 98.5
Q 76.7 93.8, 79.0 91.3
Q 81.3 88.7, 85.5 88.7
Q 89.8 88.7, 92.1 91.3
Q 94.4 93.8, 94.4 98.5
Q 94.4 103.2, 92.1 105.8
Q 89.7 108.5, 85.5 108.5
Q 81.3 108.5, 79.0 105.8
Q 76.7 103.2, 76.7 98.5
M 85.5 106.3
Q 88.5 106.3, 90.0 104.3
Q 91.7 102.3, 91.7 98.5
Q 91.7 94.7, 90.0 92.8
Q 88.5 90.9, 85.5 90.9
Q 82.6 90.9, 81.0 92.8
Q 79.4 94.7, 79.4 98.5
Q 79.4 102.4, 81.0 104.3
Q 82.6 106.3, 85.5 106.3
' fill='#FF0000'/>
<path class='atom-2' d='M 8.0 54.1
Q 8.5 52.9, 9.6 52.2
Q 10.7 51.5, 12.3 51.5
Q 14.2 51.5, 15.3 52.6
Q 16.4 53.6, 16.4 55.5
Q 16.4 57.4, 15.0 59.1
Q 13.6 60.9, 10.7 63.0
L 16.6 63.0
L 16.6 64.4
L 8.0 64.4
L 8.0 63.2
Q 10.4 61.5, 11.8 60.2
Q 13.2 59.0, 13.9 57.9
Q 14.6 56.7, 14.6 55.6
Q 14.6 54.3, 13.9 53.7
Q 13.3 53.0, 12.3 53.0
Q 11.3 53.0, 10.6 53.4
Q 9.9 53.8, 9.4 54.7
L 8.0 54.1
' fill='#000000'/>
<path class='atom-2' d='M 19.1 54.8
L 21.7 54.8
L 21.7 63.0
L 31.5 63.0
L 31.5 54.8
L 34.1 54.8
L 34.1 74.0
L 31.5 74.0
L 31.5 65.1
L 21.7 65.1
L 21.7 74.0
L 19.1 74.0
L 19.1 54.8
' fill='#000000'/>
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
<path class='bond-0 atom-0 atom-1' d='M 65.2,98.4 L 64.9,97.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 62.3,100.4 L 61.8,99.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 59.4,102.4 L 58.6,101.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 56.5,104.3 L 55.4,102.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 53.6,106.3 L 52.3,104.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 50.7,108.3 L 49.1,105.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 47.8,110.3 L 46.0,107.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 44.9,112.2 L 42.8,108.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 42.0,114.2 L 39.6,110.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 39.1,116.2 L 36.5,111.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 36.2,118.2 L 33.3,113.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 33.3,120.2 L 30.1,114.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 30.4,122.1 L 27.0,116.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 27.5,124.1 L 23.8,117.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 24.6,126.1 L 20.7,119.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 21.7,128.1 L 17.5,120.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 18.8,130.0 L 14.3,122.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 15.9,132.0 L 11.2,123.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 13.0,134.0 L 8.0,125.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:1.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 68.1,96.4 L 68.1,69.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 68.1,69.0 L 68.1,41.5' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 68.1,96.4 L 91.5,109.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 91.5,109.9 L 114.8,123.4' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 68.1,95.1 L 68.1,96.4 L 69.3,97.1' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-2' d='M 62.5 20.5
L 73.7 20.5
L 73.7 22.7
L 65.0 22.7
L 65.0 28.4
L 72.7 28.4
L 72.7 30.6
L 65.0 30.6
L 65.0 39.3
L 62.5 39.3
L 62.5 20.5
' fill='#33CCCC'/>
<path class='atom-3' d='M 117.0 129.7
Q 117.0 125.2, 119.3 122.7
Q 121.5 120.2, 125.7 120.2
Q 129.9 120.2, 132.1 122.7
Q 134.3 125.2, 134.3 129.7
Q 134.3 134.3, 132.1 136.9
Q 129.8 139.5, 125.7 139.5
Q 121.5 139.5, 119.3 136.9
Q 117.0 134.3, 117.0 129.7
M 125.7 137.4
Q 128.6 137.4, 130.1 135.5
Q 131.7 133.5, 131.7 129.7
Q 131.7 126.0, 130.1 124.2
Q 128.6 122.3, 125.7 122.3
Q 122.8 122.3, 121.2 124.1
Q 119.7 126.0, 119.7 129.7
Q 119.7 133.5, 121.2 135.5
Q 122.8 137.4, 125.7 137.4
' fill='#FF0000'/>
<path class='atom-3' d='M 137.3 120.4
L 139.8 120.4
L 139.8 128.4
L 149.4 128.4
L 149.4 120.4
L 152.0 120.4
L 152.0 139.2
L 149.4 139.2
L 149.4 130.5
L 139.8 130.5
L 139.8 139.2
L 137.3 139.2
L 137.3 120.4
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
<path class='bond-0 atom-0 atom-1' d='M 68.1,96.4 L 13.0,134.0 L 8.0,125.4 Z' style='fill:#000000;fill-rule:evenodd;fill-opacity:1;stroke:#000000;stroke-width:0.5px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='bond-1 atom-0 atom-2' d='M 68.1,96.4 L 68.1,69.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 68.1,69.0 L 68.1,41.5' style='fill:none;fill-rule:evenodd;stroke:#33CCCC;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 68.1,96.4 L 91.5,109.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 91.5,109.9 L 114.8,123.4' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 68.1,95.1 L 68.1,96.4 L 69.3,97.1' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-2' d='M 62.5 20.5
L 73.7 20.5
L 73.7 22.7
L 65.0 22.7
L 65.0 28.4
L 72.7 28.4
L 72.7 30.6
L 65.0 30.6
L 65.0 39.3
L 62.5 39.3
L 62.5 20.5
' fill='#33CCCC'/>
<path class='atom-3' d='M 117.0 129.7
Q 117.0 125.2, 119.3 122.7
Q 121.5 120.2, 125.7 120.2
Q 129.9 120.2, 132.1 122.7
Q 134.3 125.2, 134.3 129.7
Q 134.3 134.3, 132.1 136.9
Q 129.8 139.5, 125.7 139.5
Q 121.5 139.5, 119.3 136.9
Q 117.0 134.3, 117.0 129.7
M 125.7 137.4
Q 128.6 137.4, 130.1 135.5
Q 131.7 133.5, 131.7 129.7
Q 131.7 126.0, 130.1 124.2
Q 128.6 122.3, 125.7 122.3
Q 122.8 122.3, 121.2 124.1
Q 119.7 126.0, 119.7 129.7
Q 119.7 133.5, 121.2 135.5
Q 122.8 137.4, 125.7 137.4
' fill='#FF0000'/>
<path class='atom-3' d='M 137.3 120.4
L 139.8 120.4
L 139.8 128.4
L 149.4 128.4
L 149.4 120.4
L 152.0 120.4
L 152.0 139.2
L 149.4 139.2
L 149.4 130.5
L 139.8 130.5
L 139.8 139.2
L 137.3 139.2
L 137.3 120.4
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
<path class='bond-0 atom-0 atom-1' d='M 62.2,92.0 L 44.0,102.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 44.0,102.5 L 25.8,113.1' style='fill:none;fill-rule:evenodd;stroke:#7F4C19;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 62.2,92.0 L 62.2,72.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 62.2,72.9 L 62.2,53.8' style='fill:none;fill-rule:evenodd;stroke:#00CC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 62.2,92.0 L 103.0,115.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 69.3,87.9 L 103.0,107.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 103.0,115.6 L 120.0,105.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 120.0,105.7 L 137.1,95.9' style='fill:none;fill-rule:evenodd;stroke:#00CC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 101.0,114.4 L 103.0,115.6 L 103.8,115.1' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-1' d='M 14.5 115.2
Q 15.8 115.6, 16.4 116.4
Q 17.1 117.1, 17.1 118.3
Q 17.1 120.1, 15.9 121.2
Q 14.7 122.2, 12.5 122.2
L 8.0 122.2
L 8.0 108.9
L 11.9 108.9
Q 14.2 108.9, 15.4 109.8
Q 16.5 110.7, 16.5 112.4
Q 16.5 114.4, 14.5 115.2
M 9.8 110.4
L 9.8 114.6
L 11.9 114.6
Q 13.3 114.6, 13.9 114.0
Q 14.6 113.5, 14.6 112.4
Q 14.6 110.4, 11.9 110.4
L 9.8 110.4
M 12.5 120.7
Q 13.8 120.7, 14.5 120.1
Q 15.2 119.5, 15.2 118.3
Q 15.2 117.2, 14.4 116.6
Q 13.7 116.1, 12.2 116.1
L 9.8 116.1
L 9.8 120.7
L 12.5 120.7
' fill='#7F4C19'/>
<path class='atom-1' d='M 20.1 112.5
L 20.3 113.9
Q 21.3 112.4, 23.0 112.4
Q 23.5 112.4, 24.2 112.6
L 23.9 114.1
Q 23.1 114.0, 22.7 114.0
Q 21.9 114.0, 21.4 114.3
Q 20.8 114.6, 20.4 115.3
L 20.4 122.2
L 18.7 122.2
L 18.7 112.5
L 20.1 112.5
' fill='#7F4C19'/>
<path class='atom-2' d='M 57.1 45.4
Q 57.1 42.1, 58.6 40.3
Q 60.2 38.6, 63.1 38.6
Q 65.9 38.6, 67.4 40.5
L 66.1 41.5
Q 65.0 40.1, 63.1 40.1
Q 61.1 40.1, 60.0 41.5
Q 59.0 42.8, 59.0 45.4
Q 59.0 48.0, 60.1 49.4
Q 61.2 50.7, 63.3 50.7
Q 64.8 50.7, 66.5 49.9
L 67.1 51.3
Q 66.4 51.7, 65.3 52.0
Q 64.2 52.2, 63.1 52.2
Q 60.2 52.2, 58.6 50.5
Q 57.1 48.7, 57.1 45.4
' fill='#00CC00'/>
<path class='atom-2' d='M 68.7 37.8
L 70.4 37.8
L 70.4 52.1
L 68.7 52.1
L 68.7 37.8
' fill='#00CC00'/>
<path class='atom-4' d='M 138.6 92.5
Q 138.6 89.2, 140.2 87.4
Q 141.7 85.7, 144.7 85.7
Q 147.4 85.7, 148.9 87.6
L 147.7 88.6
Q 146.6 87.2, 144.7 87.2
Q 142.7 87.2, 141.6 88.6
Q 140.6 89.9, 140.6 92.5
Q 140.6 95.1, 141.6 96.5
Q 142.8 97.8, 144.9 97.8
Q 146.4 97.8, 148.1 96.9
L 148.6 98.3
Q 147.9 98.8, 146.9 99.1
Q 145.8 99.3, 144.6 99.3
Q 141.7 99.3, 140.2 97.6
Q 138.6 95.8, 138.6 92.5
' fill='#00CC00'/>
<path class='atom-4' d='M 150.3 84.9
L 152.0 84.9
L 152.0 99.2
L 150.3 99.2
L 150.3 84.9
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
<path class='bond-0 atom-0 atom-1' d='M 102.7,91.5 L 102.7,72.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 102.7,72.0 L 102.7,52.4' style='fill:none;fill-rule:evenodd;stroke:#7F4C19;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 102.7,91.5 L 119.8,101.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-0 atom-2' d='M 119.8,101.4 L 137.0,111.3' style='fill:none;fill-rule:evenodd;stroke:#00CC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 102.7,91.5 L 61.6,115.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-0 atom-3' d='M 95.6,87.4 L 61.6,107.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 61.6,115.2 L 42.3,104.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 42.3,104.1 L 23.0,92.9' style='fill:none;fill-rule:evenodd;stroke:#00CC00;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 63.7,114.0 L 61.6,115.2 L 60.7,114.7' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-1' d='M 104.6 43.8
Q 105.9 44.1, 106.6 44.9
Q 107.2 45.7, 107.2 46.9
Q 107.2 48.7, 106.0 49.8
Q 104.9 50.8, 102.6 50.8
L 98.1 50.8
L 98.1 37.4
L 102.1 37.4
Q 104.4 37.4, 105.5 38.3
Q 106.7 39.3, 106.7 41.0
Q 106.7 43.0, 104.6 43.8
M 99.9 38.9
L 99.9 43.1
L 102.1 43.1
Q 103.4 43.1, 104.1 42.6
Q 104.8 42.1, 104.8 41.0
Q 104.8 38.9, 102.1 38.9
L 99.9 38.9
M 102.6 49.3
Q 103.9 49.3, 104.6 48.7
Q 105.3 48.1, 105.3 46.9
Q 105.3 45.8, 104.6 45.2
Q 103.8 44.7, 102.3 44.7
L 99.9 44.7
L 99.9 49.3
L 102.6 49.3
' fill='#7F4C19'/>
<path class='atom-1' d='M 110.3 41.1
L 110.5 42.4
Q 111.5 40.9, 113.2 40.9
Q 113.7 40.9, 114.4 41.1
L 114.1 42.7
Q 113.3 42.5, 112.9 42.5
Q 112.1 42.5, 111.6 42.8
Q 111.0 43.1, 110.6 43.9
L 110.6 50.8
L 108.8 50.8
L 108.8 41.1
L 110.3 41.1
' fill='#7F4C19'/>
<path class='atom-2' d='M 138.5 115.7
Q 138.5 112.3, 140.1 110.6
Q 141.7 108.8, 144.6 108.8
Q 147.4 108.8, 148.9 110.8
L 147.6 111.8
Q 146.6 110.4, 144.6 110.4
Q 142.6 110.4, 141.5 111.8
Q 140.5 113.1, 140.5 115.7
Q 140.5 118.3, 141.6 119.7
Q 142.7 121.1, 144.9 121.1
Q 146.3 121.1, 148.1 120.2
L 148.6 121.6
Q 147.9 122.1, 146.8 122.3
Q 145.8 122.6, 144.6 122.6
Q 141.7 122.6, 140.1 120.8
Q 138.5 119.0, 138.5 115.7
' fill='#00CC00'/>
<path class='atom-2' d='M 150.3 108.0
L 152.0 108.0
L 152.0 122.4
L 150.3 122.4
L 150.3 108.0
' fill='#00CC00'/>
<path class='atom-4' d='M 8.0 92.0
Q 8.0 88.6, 9.6 86.9
Q 11.1 85.1, 14.1 85.1
Q 16.9 85.1, 18.4 87.1
L 17.1 88.1
Q 16.0 86.7, 14.1 86.7
Q 12.1 86.7, 11.0 88.1
Q 9.9 89.4, 9.9 92.0
Q 9.9 94.6, 11.0 96.0
Q 12.2 97.4, 14.3 97.4
Q 15.8 97.4, 17.5 96.5
L 18.0 97.9
Q 17.3 98.4, 16.3 98.6
Q 15.2 98.9, 14.0 98.9
Q 11.1 98.9, 9.6 97.1
Q 8.0 95.3, 8.0 92.0
' fill='#00CC00'/>
<path class='atom-4' d='M 19.7 84.3
L 21.5 84.3
L 21.5 98.7
L 19.7 98.7
L 19.7 84.3
' fill='#00CC00'/>
</svg>
</span>
</div>


### Groups
The following abbreviations may be used to represent various functional groups:
```py
(5aN), (5aNbN), (5aNbO), (5aNbS), (5aNcN), (5aNcO), (5aNcS), (5aNdN), (5aNdO), (5aNdS), (5aNeN), (5aNeO), (5aNeS), (5aO), (5aS), (5bN), (5bNaN), (5bNaO), (5bNaS), (5bNcN), (5bNcO), (5bNcS), (5bNdN), (5bNdO), (5bNdS), (5bNeN), (5bNeO), (5bNeS), (5bO), (5bS), (5cN), (5cNaN), (5cNaO), (5cNaS), (5cNbN), (5cNbO), (5cNbS), (5cNdN), (5cNdO), (5cNdS), (5cNeN), (5cNeO), (5cNeS), (5cO), (5cS), (5dN), (5dNaN), (5dNaO), (5dNaS), (5dNbN), (5dNbO), (5dNbS), (5dNcN), (5dNcO), (5dNcS), (5dNeN), (5dNeO), (5dNeS), (5dO), (5dS), (5eN), (5eNaN), (5eNaO), (5eNaS), (5eNbN), (5eNbO), (5eNbS), (5eNcN), (5eNcO), (5eNcS), (5eNdN), (5eNdO), (5eNdS), (5eO), (5eS), (6), (6aN), (6abN), (6abcN), (6abdN), (6abeN), (6abfN), (6acN), (6acdN), (6aceN), (6acfN), (6adN), (6adeN), (6adfN), (6aeN), (6aefN), (6afN), (6bN), (6bcN), (6bcdN), (6bceN), (6bcfN), (6bdN), (6bdeN), (6bdfN), (6beN), (6befN), (6bfN), (6cN), (6cdN), (6cdeN), (6cdfN), (6ceN), (6cefN), (6cfN), (6dN), (6deN), (6defN), (6dfN), (6eN), (6efN), (6fN), [Ac], [Bn], [Boc], [Bz], [CCl3], [CF3], [CHO], [CN], [COO-], [COOEt], [COOH], [COOMe], [Cbz], [Cy], [Et], [Ms], [NC], [NHAc], [NHMe], [NMe2], [NO2], [OAc], [OEt], [OMe], [OiBu], [PO3], [Ph], [Piv], [SMe], [SO3], [Tf], [Tol], [Ts], [iBu], [iPr], [nBu], [nDec], [nHept], [nHex], [nNon], [nOct], [nPent], [nPr], [sBu], [tBu]
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
<path class='bond-0 atom-1 atom-0' d='M 56.1,76.4 L 57.7,66.1 L 59.7,66.6 Z' style='fill:#000000;fill-rule:evenodd;fill-opacity:1;stroke:#000000;stroke-width:0.5px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='bond-0 atom-1 atom-0' d='M 57.7,66.1 L 63.3,56.7 L 59.4,55.7 Z' style='fill:#0000FF;fill-rule:evenodd;fill-opacity:1;stroke:#0000FF;stroke-width:0.5px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='bond-0 atom-1 atom-0' d='M 57.7,66.1 L 59.7,66.6 L 63.3,56.7 Z' style='fill:#0000FF;fill-rule:evenodd;fill-opacity:1;stroke:#0000FF;stroke-width:0.5px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='bond-1 atom-1 atom-2' d='M 56.1,76.4 L 75.0,95.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 75.0,95.1 L 100.7,88.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 100.7,88.1 L 107.4,62.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 105.1,86.9 L 110.7,65.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 107.4,62.3 L 133.1,55.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-6' d='M 133.1,55.3 L 152.0,74.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-6' d='M 131.9,59.7 L 147.6,75.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-6 atom-7' d='M 152.0,74.0 L 145.3,99.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-7 atom-8' d='M 145.3,99.7 L 119.6,106.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-7 atom-8' d='M 142.0,96.5 L 120.8,102.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-8 atom-3' d='M 119.6,106.8 L 100.7,88.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-1 atom-9' d='M 56.1,76.4 L 30.4,83.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 28.8,81.8 L 25.9,92.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 25.9,92.9 L 23.0,104.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 32.6,82.9 L 29.7,93.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 29.7,93.9 L 26.8,105.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-9 atom-11' d='M 30.4,83.5 L 23.1,76.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-9 atom-11' d='M 23.1,76.3 L 15.8,69.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 57.0,77.4 L 56.1,76.4 L 54.8,76.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 74.1,94.2 L 75.0,95.1 L 76.3,94.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 107.1,63.6 L 107.4,62.3 L 108.7,62.0' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 131.8,55.6 L 133.1,55.3 L 134.0,56.2' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 151.1,73.0 L 152.0,74.0 L 151.7,75.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 145.6,98.4 L 145.3,99.7 L 144.0,100.1' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 120.9,106.4 L 119.6,106.8 L 118.7,105.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 31.7,83.1 L 30.4,83.5 L 30.0,83.1' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-0' d='M 61.1 46.9
L 63.6 50.9
Q 63.8 51.3, 64.2 52.0
Q 64.6 52.7, 64.7 52.8
L 64.7 46.9
L 65.7 46.9
L 65.7 54.4
L 64.6 54.4
L 62.0 50.1
Q 61.7 49.6, 61.3 49.0
Q 61.0 48.4, 60.9 48.2
L 60.9 54.4
L 59.9 54.4
L 59.9 46.9
L 61.1 46.9
' fill='#0000FF'/>
<path class='atom-0' d='M 67.1 46.9
L 68.1 46.9
L 68.1 50.1
L 72.0 50.1
L 72.0 46.9
L 73.0 46.9
L 73.0 54.4
L 72.0 54.4
L 72.0 51.0
L 68.1 51.0
L 68.1 54.4
L 67.1 54.4
L 67.1 46.9
' fill='#0000FF'/>
<path class='atom-0' d='M 76.8 55.6
Q 77.3 55.7, 77.5 56.0
Q 77.8 56.3, 77.8 56.8
Q 77.8 57.3, 77.6 57.6
Q 77.3 57.9, 77.0 58.1
Q 76.6 58.3, 76.1 58.3
Q 75.5 58.3, 75.1 58.1
Q 74.7 57.9, 74.4 57.5
L 74.8 57.1
Q 75.1 57.5, 75.4 57.6
Q 75.7 57.7, 76.1 57.7
Q 76.5 57.7, 76.8 57.5
Q 77.1 57.2, 77.1 56.8
Q 77.1 56.3, 76.8 56.1
Q 76.5 55.9, 75.9 55.9
L 75.5 55.9
L 75.5 55.4
L 75.8 55.4
Q 76.4 55.4, 76.7 55.1
Q 77.0 54.9, 77.0 54.4
Q 77.0 54.1, 76.7 53.9
Q 76.5 53.7, 76.1 53.7
Q 75.7 53.7, 75.4 53.9
Q 75.1 54.0, 74.9 54.4
L 74.5 54.1
Q 74.6 53.7, 75.1 53.5
Q 75.5 53.2, 76.1 53.2
Q 76.8 53.2, 77.2 53.5
Q 77.7 53.9, 77.7 54.4
Q 77.7 54.8, 77.4 55.1
Q 77.2 55.4, 76.8 55.6
' fill='#0000FF'/>
<path class='atom-0' d='M 74.5 48.3
L 75.8 48.3
L 75.8 46.9
L 76.4 46.9
L 76.4 48.3
L 77.7 48.3
L 77.7 48.8
L 76.4 48.8
L 76.4 50.2
L 75.8 50.2
L 75.8 48.8
L 74.5 48.8
L 74.5 48.3
' fill='#0000FF'/>
<path class='atom-10' d='M 20.2 109.2
Q 20.2 107.4, 21.1 106.4
Q 22.0 105.4, 23.7 105.4
Q 25.3 105.4, 26.2 106.4
Q 27.1 107.4, 27.1 109.2
Q 27.1 111.1, 26.2 112.1
Q 25.3 113.1, 23.7 113.1
Q 22.0 113.1, 21.1 112.1
Q 20.2 111.1, 20.2 109.2
M 23.7 112.3
Q 24.8 112.3, 25.4 111.5
Q 26.1 110.7, 26.1 109.2
Q 26.1 107.8, 25.4 107.0
Q 24.8 106.3, 23.7 106.3
Q 22.5 106.3, 21.9 107.0
Q 21.3 107.7, 21.3 109.2
Q 21.3 110.8, 21.9 111.5
Q 22.5 112.3, 23.7 112.3
' fill='#FF0000'/>
<path class='atom-11' d='M 8.0 64.8
Q 8.0 63.0, 8.9 62.0
Q 9.8 60.9, 11.5 60.9
Q 13.1 60.9, 14.0 62.0
Q 14.9 63.0, 14.9 64.8
Q 14.9 66.6, 14.0 67.7
Q 13.1 68.7, 11.5 68.7
Q 9.8 68.7, 8.9 67.7
Q 8.0 66.6, 8.0 64.8
M 11.5 67.8
Q 12.6 67.8, 13.2 67.1
Q 13.9 66.3, 13.9 64.8
Q 13.9 63.3, 13.2 62.6
Q 12.6 61.8, 11.5 61.8
Q 10.3 61.8, 9.7 62.5
Q 9.1 63.3, 9.1 64.8
Q 9.1 66.3, 9.7 67.1
Q 10.3 67.8, 11.5 67.8
' fill='#FF0000'/>
<path class='atom-11' d='M 16.0 62.2
L 18.6 62.2
L 18.6 62.7
L 16.0 62.7
L 16.0 62.2
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
<path class='bond-0 atom-0 atom-1' d='M 91.4,85.7 L 71.2,86.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 89.5,82.8 L 72.8,83.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 71.2,86.9 L 60.0,70.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 60.0,70.1 L 69.1,52.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 63.5,69.9 L 71.0,54.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 69.1,52.0 L 89.3,50.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 89.3,50.7 L 100.4,67.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 87.7,53.9 L 96.9,67.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-0' d='M 100.4,67.6 L 91.4,85.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-5 atom-6' d='M 100.4,67.6 L 120.6,66.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-6 atom-7' d='M 120.6,66.4 L 131.8,83.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-7 atom-8' d='M 131.8,83.2 L 152.0,82.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-7 atom-9' d='M 131.8,83.2 L 122.8,101.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-2 atom-10' d='M 60.0,70.1 L 39.9,71.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-10 atom-11' d='M 39.9,71.3 L 28.7,54.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-12 atom-10 atom-12' d='M 39.9,71.3 L 30.8,89.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-13 atom-12 atom-13' d='M 31.6,87.9 L 22.7,88.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-13 atom-12 atom-13' d='M 22.7,88.4 L 13.8,88.9' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-13 atom-12 atom-13' d='M 31.8,90.9 L 22.9,91.4' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-13 atom-12 atom-13' d='M 22.9,91.4 L 14.0,92.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-14 atom-12 atom-14' d='M 30.8,89.4 L 35.2,96.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-14 atom-12 atom-14' d='M 35.2,96.1 L 39.6,102.7' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 90.4,85.8 L 91.4,85.7 L 91.9,84.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 72.2,86.9 L 71.2,86.9 L 70.7,86.1' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 68.6,52.9 L 69.1,52.0 L 70.1,51.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 88.3,50.8 L 89.3,50.7 L 89.8,51.6' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 119.6,66.4 L 120.6,66.4 L 121.2,67.2' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 31.3,88.5 L 30.8,89.4 L 31.0,89.8' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-13' d='M 8.0 90.7
Q 8.0 89.3, 8.7 88.5
Q 9.4 87.8, 10.6 87.8
Q 11.9 87.8, 12.6 88.5
Q 13.3 89.3, 13.3 90.7
Q 13.3 92.1, 12.6 92.9
Q 11.9 93.6, 10.6 93.6
Q 9.4 93.6, 8.7 92.9
Q 8.0 92.1, 8.0 90.7
M 10.6 93.0
Q 11.5 93.0, 12.0 92.4
Q 12.5 91.8, 12.5 90.7
Q 12.5 89.6, 12.0 89.0
Q 11.5 88.4, 10.6 88.4
Q 9.8 88.4, 9.3 89.0
Q 8.8 89.5, 8.8 90.7
Q 8.8 91.8, 9.3 92.4
Q 9.8 93.0, 10.6 93.0
' fill='#FF0000'/>
<path class='atom-14' d='M 39.4 106.3
Q 39.4 104.9, 40.0 104.2
Q 40.7 103.4, 42.0 103.4
Q 43.3 103.4, 43.9 104.2
Q 44.6 104.9, 44.6 106.3
Q 44.6 107.7, 43.9 108.5
Q 43.3 109.3, 42.0 109.3
Q 40.7 109.3, 40.0 108.5
Q 39.4 107.7, 39.4 106.3
M 42.0 108.6
Q 42.9 108.6, 43.3 108.0
Q 43.8 107.5, 43.8 106.3
Q 43.8 105.2, 43.3 104.6
Q 42.9 104.0, 42.0 104.0
Q 41.1 104.0, 40.6 104.6
Q 40.2 105.2, 40.2 106.3
Q 40.2 107.5, 40.6 108.0
Q 41.1 108.6, 42.0 108.6
' fill='#FF0000'/>
<path class='atom-14' d='M 45.5 103.5
L 46.3 103.5
L 46.3 105.9
L 49.2 105.9
L 49.2 103.5
L 50.0 103.5
L 50.0 109.2
L 49.2 109.2
L 49.2 106.5
L 46.3 106.5
L 46.3 109.2
L 45.5 109.2
L 45.5 103.5
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
<path class='bond-0 atom-0 atom-1' d='M 77.7,82.8 L 102.5,77.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 102.5,77.8 L 119.2,96.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 106.8,77.0 L 120.6,92.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-2 atom-3' d='M 119.2,96.8 L 143.9,91.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 143.9,91.8 L 152.0,67.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 141.0,88.5 L 147.7,68.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 152.0,67.9 L 135.3,48.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-6' d='M 135.3,48.9 L 110.6,53.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-5 atom-6' d='M 133.9,53.1 L 113.5,57.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-6 atom-1' d='M 110.6,53.9 L 102.5,77.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-0 atom-7' d='M 77.7,82.8 L 86.9,106.3' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-7 atom-8' d='M 86.9,106.3 L 71.1,126.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-0 atom-9' d='M 77.7,82.8 L 77.2,57.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 77.2,59.7 L 86.4,54.1' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 86.4,54.1 L 95.6,48.5' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 75.2,56.5 L 84.4,50.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-10 atom-9 atom-10' d='M 84.4,50.9 L 93.6,45.3' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-9 atom-11' d='M 77.2,57.5 L 67.8,52.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-11 atom-9 atom-11' d='M 67.8,52.5 L 58.5,47.4' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-12 atom-11 atom-12' d='M 51.5,47.6 L 42.4,53.1' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-12 atom-11 atom-12' d='M 42.4,53.1 L 33.4,58.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-13 atom-12 atom-13' d='M 35.3,57.4 L 25.8,52.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-13 atom-12 atom-13' d='M 25.8,52.2 L 16.3,47.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-13 atom-12 atom-13' d='M 33.5,60.7 L 24.0,55.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-13 atom-12 atom-13' d='M 24.0,55.5 L 14.5,50.4' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-14 atom-12 atom-14' d='M 33.4,58.6 L 33.7,69.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-14 atom-12 atom-14' d='M 33.7,69.0 L 33.9,79.4' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-15 atom-14 atom-15' d='M 37.6,85.7 L 46.9,90.8' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-15 atom-14 atom-15' d='M 46.9,90.8 L 56.2,95.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-16 atom-15 atom-16' d='M 54.3,94.9 L 54.5,105.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-16 atom-15 atom-16' d='M 54.5,105.8 L 54.8,116.7' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-16 atom-15 atom-16' d='M 58.0,94.8 L 58.3,105.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-16 atom-15 atom-16' d='M 58.3,105.7 L 58.5,116.7' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-17 atom-15 atom-0' d='M 56.2,95.9 L 77.7,82.8' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 118.3,95.8 L 119.2,96.8 L 120.4,96.5' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 142.7,92.1 L 143.9,91.8 L 144.3,90.6' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 151.6,69.1 L 152.0,67.9 L 151.2,67.0' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 136.2,49.9 L 135.3,48.9 L 134.1,49.2' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 111.8,53.7 L 110.6,53.9 L 110.2,55.1' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 86.5,105.1 L 86.9,106.3 L 86.1,107.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 77.2,58.8 L 77.2,57.5 L 76.7,57.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 33.9,58.3 L 33.4,58.6 L 33.4,59.1' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 55.7,95.7 L 56.2,95.9 L 57.2,95.3' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-10' d='M 95.4 44.4
Q 95.4 42.7, 96.3 41.8
Q 97.1 40.8, 98.7 40.8
Q 100.3 40.8, 101.1 41.8
Q 102.0 42.7, 102.0 44.4
Q 102.0 46.2, 101.1 47.2
Q 100.3 48.1, 98.7 48.1
Q 97.1 48.1, 96.3 47.2
Q 95.4 46.2, 95.4 44.4
M 98.7 47.3
Q 99.8 47.3, 100.4 46.6
Q 101.0 45.9, 101.0 44.4
Q 101.0 43.0, 100.4 42.3
Q 99.8 41.6, 98.7 41.6
Q 97.6 41.6, 97.0 42.3
Q 96.4 43.0, 96.4 44.4
Q 96.4 45.9, 97.0 46.6
Q 97.6 47.3, 98.7 47.3
' fill='#FF0000'/>
<path class='atom-11' d='M 53.4 41.9
L 55.8 45.6
Q 56.0 46.0, 56.4 46.7
Q 56.7 47.4, 56.8 47.4
L 56.8 41.9
L 57.7 41.9
L 57.7 49.0
L 56.7 49.0
L 54.2 44.9
Q 53.9 44.4, 53.6 43.8
Q 53.3 43.3, 53.2 43.1
L 53.2 49.0
L 52.3 49.0
L 52.3 41.9
L 53.4 41.9
' fill='#0000FF'/>
<path class='atom-11' d='M 52.2 34.0
L 53.2 34.0
L 53.2 37.0
L 56.8 37.0
L 56.8 34.0
L 57.8 34.0
L 57.8 41.1
L 56.8 41.1
L 56.8 37.8
L 53.2 37.8
L 53.2 41.1
L 52.2 41.1
L 52.2 34.0
' fill='#0000FF'/>
<path class='atom-13' d='M 8.0 46.5
Q 8.0 44.7, 8.8 43.8
Q 9.7 42.8, 11.3 42.8
Q 12.9 42.8, 13.7 43.8
Q 14.6 44.7, 14.6 46.5
Q 14.6 48.2, 13.7 49.2
Q 12.8 50.2, 11.3 50.2
Q 9.7 50.2, 8.8 49.2
Q 8.0 48.2, 8.0 46.5
M 11.3 49.4
Q 12.4 49.4, 13.0 48.6
Q 13.6 47.9, 13.6 46.5
Q 13.6 45.1, 13.0 44.4
Q 12.4 43.6, 11.3 43.6
Q 10.2 43.6, 9.6 44.3
Q 9.0 45.0, 9.0 46.5
Q 9.0 47.9, 9.6 48.6
Q 10.2 49.4, 11.3 49.4
' fill='#FF0000'/>
<path class='atom-14' d='M 24.3 80.2
L 25.3 80.2
L 25.3 83.3
L 28.9 83.3
L 28.9 80.2
L 29.9 80.2
L 29.9 87.4
L 28.9 87.4
L 28.9 84.1
L 25.3 84.1
L 25.3 87.4
L 24.3 87.4
L 24.3 80.2
' fill='#0000FF'/>
<path class='atom-14' d='M 32.4 80.2
L 34.8 84.0
Q 35.0 84.4, 35.4 85.1
Q 35.8 85.7, 35.8 85.8
L 35.8 80.2
L 36.7 80.2
L 36.7 87.4
L 35.8 87.4
L 33.2 83.2
Q 32.9 82.7, 32.6 82.2
Q 32.3 81.6, 32.2 81.5
L 32.2 87.4
L 31.3 87.4
L 31.3 80.2
L 32.4 80.2
' fill='#0000FF'/>
<path class='atom-16' d='M 53.5 121.2
Q 53.5 119.5, 54.3 118.5
Q 55.2 117.5, 56.8 117.5
Q 58.3 117.5, 59.2 118.5
Q 60.0 119.5, 60.0 121.2
Q 60.0 122.9, 59.2 123.9
Q 58.3 124.9, 56.8 124.9
Q 55.2 124.9, 54.3 123.9
Q 53.5 122.9, 53.5 121.2
M 56.8 124.1
Q 57.8 124.1, 58.4 123.3
Q 59.0 122.6, 59.0 121.2
Q 59.0 119.8, 58.4 119.1
Q 57.8 118.3, 56.8 118.3
Q 55.7 118.3, 55.1 119.0
Q 54.5 119.8, 54.5 121.2
Q 54.5 122.6, 55.1 123.3
Q 55.7 124.1, 56.8 124.1
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
<path class='bond-0 atom-0 atom-1' d='M 43.4,140.2 L 43.4,124.4' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 43.4,124.4 L 43.4,108.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 48.9,140.2 L 48.9,124.4' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-0 atom-0 atom-1' d='M 48.9,124.4 L 48.9,108.7' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-1 atom-1 atom-2' d='M 46.1,110.2 L 14.6,92.0' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-1 atom-3' d='M 46.1,110.2 L 58.9,102.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-2 atom-1 atom-3' d='M 58.9,102.9 L 71.7,95.5' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 83.6,95.5 L 96.4,102.9' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-3 atom-3 atom-4' d='M 96.4,102.9 L 109.1,110.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 109.1,113.4 L 122.6,105.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 122.6,105.6 L 136.1,97.8' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 106.4,108.7 L 119.9,100.9' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-4 atom-4 atom-5' d='M 119.9,100.9 L 133.4,93.1' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-4 atom-6' d='M 109.1,110.2 L 109.1,125.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-5 atom-4 atom-6' d='M 109.1,125.2 L 109.1,140.2' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-6 atom-7 atom-8' d='M 33.5,67.8 L 65.0,49.6' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-8 atom-9' d='M 67.8,51.2 L 67.8,35.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-8 atom-9' d='M 67.8,35.5 L 67.8,19.8' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-8 atom-9' d='M 62.3,51.2 L 62.3,35.5' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-7 atom-8 atom-9' d='M 62.3,35.5 L 62.3,19.8' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-8 atom-10' d='M 65.0,49.6 L 78.2,57.2' style='fill:none;fill-rule:evenodd;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-8 atom-8 atom-10' d='M 78.2,57.2 L 91.4,64.8' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-10 atom-11' d='M 101.7,64.8 L 111.9,58.9' style='fill:none;fill-rule:evenodd;stroke:#0000FF;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path class='bond-9 atom-10 atom-11' d='M 111.9,58.9 L 122.1,53.0' style='fill:none;fill-rule:evenodd;stroke:#FF0000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1' />
<path d='M 44.5,109.3 L 46.1,110.2 L 46.8,109.9' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 108.5,109.9 L 109.1,110.2 L 109.1,111.0' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path d='M 63.5,50.5 L 65.0,49.6 L 65.7,50.0' style='fill:none;stroke:#000000;stroke-width:2.0px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;' />
<path class='atom-0' d='M 41.4 146.7
Q 41.4 144.2, 42.6 142.8
Q 43.8 141.4, 46.1 141.4
Q 48.4 141.4, 49.6 142.8
Q 50.9 144.2, 50.9 146.7
Q 50.9 149.2, 49.6 150.6
Q 48.4 152.0, 46.1 152.0
Q 43.9 152.0, 42.6 150.6
Q 41.4 149.2, 41.4 146.7
M 46.1 150.8
Q 47.7 150.8, 48.5 149.8
Q 49.4 148.7, 49.4 146.7
Q 49.4 144.6, 48.5 143.6
Q 47.7 142.6, 46.1 142.6
Q 44.6 142.6, 43.7 143.6
Q 42.8 144.6, 42.8 146.7
Q 42.8 148.7, 43.7 149.8
Q 44.6 150.8, 46.1 150.8
' fill='#FF0000'/>
<path class='atom-3' d='M 72.9 92.1
Q 72.9 89.6, 74.1 88.2
Q 75.3 86.8, 77.6 86.8
Q 79.9 86.8, 81.1 88.2
Q 82.4 89.6, 82.4 92.1
Q 82.4 94.6, 81.1 96.0
Q 79.9 97.4, 77.6 97.4
Q 75.4 97.4, 74.1 96.0
Q 72.9 94.6, 72.9 92.1
M 77.6 96.3
Q 79.2 96.3, 80.1 95.2
Q 80.9 94.1, 80.9 92.1
Q 80.9 90.1, 80.1 89.0
Q 79.2 88.0, 77.6 88.0
Q 76.1 88.0, 75.2 89.0
Q 74.4 90.0, 74.4 92.1
Q 74.4 94.2, 75.2 95.2
Q 76.1 96.3, 77.6 96.3
' fill='#FF0000'/>
<path class='atom-5' d='M 135.9 92.1
Q 135.9 89.6, 137.2 88.2
Q 138.4 86.8, 140.7 86.8
Q 142.9 86.8, 144.2 88.2
Q 145.4 89.6, 145.4 92.1
Q 145.4 94.6, 144.2 96.0
Q 142.9 97.4, 140.7 97.4
Q 138.4 97.4, 137.2 96.0
Q 135.9 94.6, 135.9 92.1
M 140.7 96.3
Q 142.2 96.3, 143.1 95.2
Q 143.9 94.1, 143.9 92.1
Q 143.9 90.1, 143.1 89.0
Q 142.2 88.0, 140.7 88.0
Q 139.1 88.0, 138.2 89.0
Q 137.4 90.0, 137.4 92.1
Q 137.4 94.2, 138.2 95.2
Q 139.1 96.3, 140.7 96.3
' fill='#FF0000'/>
<path class='atom-6' d='M 104.4 146.7
Q 104.4 144.2, 105.6 142.8
Q 106.9 141.4, 109.1 141.4
Q 111.4 141.4, 112.7 142.8
Q 113.9 144.2, 113.9 146.7
Q 113.9 149.2, 112.6 150.6
Q 111.4 152.0, 109.1 152.0
Q 106.9 152.0, 105.6 150.6
Q 104.4 149.2, 104.4 146.7
M 109.1 150.8
Q 110.7 150.8, 111.6 149.8
Q 112.4 148.7, 112.4 146.7
Q 112.4 144.6, 111.6 143.6
Q 110.7 142.6, 109.1 142.6
Q 107.6 142.6, 106.7 143.6
Q 105.9 144.6, 105.9 146.7
Q 105.9 148.7, 106.7 149.8
Q 107.6 150.8, 109.1 150.8
' fill='#FF0000'/>
<path class='atom-6' d='M 115.5 141.5
L 116.9 141.5
L 116.9 145.9
L 122.1 145.9
L 122.1 141.5
L 123.5 141.5
L 123.5 151.8
L 122.1 151.8
L 122.1 147.1
L 116.9 147.1
L 116.9 151.8
L 115.5 151.8
L 115.5 141.5
' fill='#FF0000'/>
<path class='atom-9' d='M 60.3 13.2
Q 60.3 10.8, 61.5 9.4
Q 62.7 8.0, 65.0 8.0
Q 67.3 8.0, 68.5 9.4
Q 69.8 10.8, 69.8 13.2
Q 69.8 15.7, 68.5 17.2
Q 67.3 18.6, 65.0 18.6
Q 62.8 18.6, 61.5 17.2
Q 60.3 15.8, 60.3 13.2
M 65.0 17.4
Q 66.6 17.4, 67.4 16.4
Q 68.3 15.3, 68.3 13.2
Q 68.3 11.2, 67.4 10.2
Q 66.6 9.2, 65.0 9.2
Q 63.5 9.2, 62.6 10.2
Q 61.8 11.2, 61.8 13.2
Q 61.8 15.3, 62.6 16.4
Q 63.5 17.4, 65.0 17.4
' fill='#FF0000'/>
<path class='atom-10' d='M 94.3 62.6
L 97.6 68.1
Q 98.0 68.6, 98.5 69.6
Q 99.1 70.6, 99.1 70.6
L 99.1 62.6
L 100.5 62.6
L 100.5 72.9
L 99.0 72.9
L 95.4 67.0
Q 95.0 66.3, 94.5 65.5
Q 94.1 64.7, 94.0 64.4
L 94.0 72.9
L 92.6 72.9
L 92.6 62.6
L 94.3 62.6
' fill='#0000FF'/>
<path class='atom-10' d='M 92.5 74.0
L 93.9 74.0
L 93.9 78.4
L 99.2 78.4
L 99.2 74.0
L 100.6 74.0
L 100.6 84.3
L 99.2 84.3
L 99.2 79.5
L 93.9 79.5
L 93.9 84.3
L 92.5 84.3
L 92.5 74.0
' fill='#0000FF'/>
<path class='atom-11' d='M 123.3 49.6
Q 123.3 47.2, 124.5 45.8
Q 125.8 44.4, 128.1 44.4
Q 130.3 44.4, 131.6 45.8
Q 132.8 47.2, 132.8 49.6
Q 132.8 52.1, 131.5 53.6
Q 130.3 55.0, 128.1 55.0
Q 125.8 55.0, 124.5 53.6
Q 123.3 52.1, 123.3 49.6
M 128.1 53.8
Q 129.6 53.8, 130.5 52.8
Q 131.3 51.7, 131.3 49.6
Q 131.3 47.6, 130.5 46.6
Q 129.6 45.6, 128.1 45.6
Q 126.5 45.6, 125.6 46.6
Q 124.8 47.6, 124.8 49.6
Q 124.8 51.7, 125.6 52.8
Q 126.5 53.8, 128.1 53.8
' fill='#FF0000'/>
<path class='atom-11' d='M 134.4 44.5
L 135.8 44.5
L 135.8 48.9
L 141.1 48.9
L 141.1 44.5
L 142.5 44.5
L 142.5 54.8
L 141.1 54.8
L 141.1 50.0
L 135.8 50.0
L 135.8 54.8
L 134.4 54.8
L 134.4 44.5
' fill='#FF0000'/>
</svg>
</span>
</div>


## Developing
This repo uses pre-commit, so after cloning run `pip install -r requirements.txt` and `pre-commit install` prior to committing.
