## METASET: Exploring Shape and Property Spaces for Data-Driven Metamaterials Design

This repository contains the **data of shape and property diverse subsets of metamaterial unit cells** autonomously selected by METASET, a methodology that uses similarity metrics to jointly measure the diversity of unit cells in both shape and property space through Determinantal Point Processes and an efficient, greedy subset selection algorithm. In our papers, we demonstrated that a smaller yet diverse set of unit cells can lead to scalable search and unbiased learning. Our diverse datasets can be directly used by metamaterials designers, especially in data-driven design or, in the 3D case, design of functionally graded materials.

### About the Data
- 2D diverse 50x50 pixelated unit cells:
  - Sizes of subsets: 20, 100
  - Formats: PNG, TXT
  - See [./2D_diverse/README](https://github.com/lychan110/metaset/tree/master/2D_diverse) for more information.
- 3D diverse unit cell families:
  - Size of subsets: 10
  - Formats: PNG, STL, level-set functions and shape generation scripts, properties
  - See [./3D_diverse/README](https://github.com/lychan110/metaset/tree/master/3D_diverse) for more information.

### Citation
If our data has been useful in your research, please cite our work:

**Journal paper:**

Chan, Y.-C., Ahmed, F., Wang, L., and Chen, W., 2020. "METASET: Exploring Shape and Property Spaces for Data-Driven Metamaterials Design."  _Journal of Mechanical Design._ March 2021; 143(3): 031707.

- _Note:_ The journal version includes a new, more complex 2D example, and fixes the errata in the conference version.
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

**Conference paper:**

Chan, Y.-C., Ahmed, F., Wang, L., and Chen, W., 2020, "METASET: An automated data selection method for scalable data-driven design of metamaterials." _ASME 2020 Design Automation Conference_. Paper Number: DETC2020-22681.

- _Errata:_ The final structures in Fig. 6 are erroneously plotted upside down. In Table 1, "Bridge" refers to the cantilever example.

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
