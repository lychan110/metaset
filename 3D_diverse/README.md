### About the 3D Data

We provide two versions of our 3D data, one from our conference paper in _IDETC 2020_, and from our journal paper in _JMD_. In _JMD_, we filtered our dataset to include only families with a wider range of structurally feasible unit cells. Both were selected using the Hausdorff distance for shape similarity. Please see the sections below and our papers for more details.

- 3D diverse unit cell families:
  - Size of subsets: 10
  - Formats: PNG, STL, level-set functions
  - Options:
    - One sample (median unit cell) per family is provided here.
    - 100 samples per family are available at [this link for 2020_IDETC](https://northwestern.box.com/s/l3prbvm90ou0idp1s117mk7cj9ilxb6j). The zip archive of 3,000 STL files is ~930MB. If the link does not work, please contact weichen (at) northwestern (dot) edu. (Will be available for `2020_JMD` in the future.)
    - A PDF file containing the level-set functions (and crystallographic information) used to generate each family is also available. More samples can be generated from the level-set functions. (Will be available for `2020_JMD` in the future.)

#### Summary of Each 3D Subset

##### 2020_IDETC

**Size 10: Property Diverse** (median samples shown)
![P20](https://github.com/lychan110/metaset/blob/master/3D_diverse/2020_IDETC/images/3D_P10_summary.png)

**Size 10: Shape and Property Diverse** (median samples shown)
![SP20](https://github.com/lychan110/metaset/blob/master/3D_diverse/2020_IDETC/images/3D_SP10_summary.png)

**Size 10: Shape Diverse** (median samples shown)
![S20](https://github.com/lychan110/metaset/blob/master/3D_diverse/2020_IDETC/images/3D_S10_summary.png)

- _Errata:_ An error during the data generation stage caused some unit cells to become "off-center", i.e., not symmetric. A few PNG and STL may be incorrect, but the level-set functions are accurate. This is fixed in the `2020_JMD` data set.

Please cite:
`Chan, Y.-C., Ahmed, F., Wang, L., and Chen, W., 2020, "METASET: An automated data selection method for scalable data-driven design of metamaterials." _ASME 2020 Design Automation Conference_. Paper Number: DETC2020-22681.`

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

##### 2020_JMD

**Size 10: Property Diverse** (median samples shown)
![P20](https://github.com/lychan110/metaset/blob/master/3D_diverse/2020_JMD/images/3D_P10_summary_JMD.png)

**Size 10: Shape and Property Diverse** (median samples shown)
![SP20](https://github.com/lychan110/metaset/blob/master/3D_diverse/2020_JMD/images/3D_SP10_summary_JMD.png)

**Size 10: Shape Diverse** (median samples shown)
![S20](https://github.com/lychan110/metaset/blob/master/3D_diverse/2020_JMD/images/3D_S10_summary_JMD.png)

Please cite:
`Chan, Y.-C., Ahmed, F., Wang, L., and Chen, W., 2020. "METASET: Exploring Shape and Property Spaces for Data-Driven Metamaterials Design.“  _Journal of Mechanical Design_ March 2021; 143(3): 031707.`

- [arXiv version available](https://arxiv.org/abs/2006.02142)

```BibTeX
@article{Chan2020jmd,
  doi = {10.1115/1.4048629},
  url = {https://doi.org/10.1115/1.4048629},
  year = {2020},
  month = nov,
  publisher = {{ASME} International},
  volume = {143},
  number = {3},
  author = {Yu-Chin Chan and Faez Ahmed and Liwei Wang and Wei Chen},
  title = {{METASET}: Exploring Shape and Property Spaces for Data-Driven Metamaterials Design},
  journal = {Journal of Mechanical Design}
}
```