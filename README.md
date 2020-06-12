## METASET: Exploring Shape and Property Spaces for Data-Driven Metamaterials Design

This repository contains the **data of shape and property diverse subsets of metamaterial unit cells** autonomously selected by METASET, a methodology that uses similarity metrics to jointly measure the diversity of unit cells in both shape and property space through Determinantal Point Processes and an efficient, greedy subset selection algorithm. In our papers, we demonstrated that a smaller yet diverse set of unit cells can lead to scalable search and unbiased learning. Our diverse datasets can be directly used by metamaterials designers, especially in data-driven design or, in the 3D case, design of functionally graded materials.

### About the Data
- 2D diverse 50x50 pixelated unit cells:
  - Sizes of subsets: 20, 100
  - Formats: PNG, TXT of binary images.
  - See [./2D_diverse/README](https://github.com/lychan110/metaset/tree/master/2D_diverse) for more information.
- 3D diverse unit cell families:
  - Size of subsets: 10
  - Formats: PNG, STL, level-set functions
  - See [./3D_diverse/README](https://github.com/lychan110/metaset/tree/master/3D_diverse) for more information.

### Citation
If our data has been useful in your research, please cite our work:

**Journal preprint:**

Chan, Y.-C., Ahmed, F., Wang, L., and Chen, W., 2020. "METASET: Exploring Shape and Property Spaces for Data-Driven Metamaterials Design.“ [arXiv preprint arXiv:2006.02142](https://arxiv.org/abs/2006.02142). (Submitted to _Journal of Mechanical Design_).

- _Note:_ The journal version includes a new, more complex 2D example, and fixes the errata in the conference version.

**Conference paper:**

Chan, Y.-C., Ahmed, F., Wang, L., and Chen, W., 2020, "METASET: An automated data selection method for scalable data-driven design of metamaterials." _ASME 2020 Design Automation Conference_. Paper Number: DETC2020-22681.

- _Errata:_ The final structures in Fig. 6 are erroneously plotted upside down. In Table 1, "Bridge" refers to the cantilever example.

```BibTeX
@inproceedings{Chan2020idetc
  year = {2020},
  author = {Yu-Chin Chan and Faez Ahmed and Liwei Wang and Wei Chen},
  title = {METASET: An automated data selection method for scalable data-driven design of metamaterials},
  booktitle={ASME 2020 International Design Engineering Technical Conferences \& Computers and Information in Engineering Conference},
  month={August},
  address={St. Louis, MO, USA}
}
```
