### About the 2D Data
- 2D diverse 50x50 pixelated unit cells:
  - Subset sizes: 20, 100
  - Formats: PNG, TXT

**Instructions for reading the TXT format:** Each unit cell is represented as a 50 x 50 binary matrix, where 1 = solid pixel and 0 = void pixel. The first entry (first row, first column) corresponds to the upper left corner of the unit cell. For example, the unit cell below can be imported and plotted in MATLAB as follows:

<img src="https://github.com/lychan110/metaset/blob/master/2D_diverse/images/read_demo.png" width="200">

```MATLAB
Phi = readmatrix('S20_001.txt');
figure; imagesc(Phi); colormap(flipud(gray)); axis tight equal ij;
```

#### Summary of Each 2D Subset

**Size 20: Property Diverse**
![P20](https://github.com/lychan110/metaset/blob/master/2D_diverse/images/2D_P20_summary.png)

**Size 20: Shape and Property Diverse**
![SP20](https://github.com/lychan110/metaset/blob/master/2D_diverse/images/2D_SP20_summary.png)

**Size 20: Shape Diverse**
![S20](https://github.com/lychan110/metaset/blob/master/2D_diverse/images/2D_S20_summary.png)

**Size 100: Property Diverse**
![P100](https://github.com/lychan110/metaset/blob/master/2D_diverse/images/2D_P100_summary.png)

**Size 100: Shape and Property Diverse**
![P100](https://github.com/lychan110/metaset/blob/master/2D_diverse/images/2D_SP100_summary.png)

**Size 100: Shape Diverse**
![P100](https://github.com/lychan110/metaset/blob/master/2D_diverse/images/2D_S100_summary.png)

Please cite:
Chan, Y.-C., Ahmed, F., Wang, L., and Chen, W., 2020, "METASET: An automated data selection method for scalable data-driven design of metamaterials." _ASME 2020 Design Automation Conference_. Paper Number: DETC2020-22681.

```BibTeX
@inproceedings{Chan2020idetc,
  doi = {10.1115/detc2020-22681},
  url = {https://doi.org/10.1115/detc2020-22681},
  year = {2020},
  month = aug,
  publisher = {American Society of Mechanical Engineers},
  author = {Yu-Chin Chan and Faez Ahmed and Liwei Wang and Wei Chen},
  title = {{METASET}: An Automated Data Selection Method for Scalable Data-Driven Design of Metamaterials},
  booktitle = {Volume 11A: 46th Design Automation Conference ({DAC})},
  address={St. Louis, MO, USA}
}
```