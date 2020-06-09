## METASET: Exploring Shape and Property Spaces for Data-Driven Metamaterials Design

This repository contains the **data of shape and property diverse subsets of metamaterial unit cells** autonomously selected by METASET, a methodology that uses similarity metrics to jointly measure the diversity of unit cells in both shape and property space through Determinantal Point Processes and an efficient, greedy subset selection algorithm. In our papers, we demonstrated that a smaller yet diverse set of unit cells can lead to scalable search and unbiased learning. Our diverse datasets can be directly used by metamaterials designers, especially in data-driven design or, in the 3D case, design of functionally graded materials.

### About the Data
- 2D diverse 50x50 pixelated unit cells:
  - Subset sizes: 20, 100
  - Formats: PNG, TXT of binary images, where 1 = solid pixel and 0 = void pixel
- 3D diverse unit cell families:
  - Subset size: 10
  - Formats: PNG, STL
  - Note: one sample (median unit cell) per family is provided. More samples can be generated from the level-set functions.
  - Additional information: PDF file containing the level-set functions used to generate each family (COMING SOON)

### Citation
If our data has been useful in your research, please do cite our work:

#### Journal preprint:
_Note:_ The journal version includes an additional 2D example, and fixes the errors in the conference version.

Chan, Y.-C., Ahmed, F., Wang, L., and Chen, W., 2020. "METASET: Exploring Shape and Property Spaces for Data-Driven Metamaterials Design.“ [arXiv preprint arXiv:2006.02142](https://arxiv.org/abs/2006.02142). (Submitted to _Journal of Mechanical Design_).

#### Conference paper:
_Errata:_ The final structures in Fig. 6 are erroneously plotted upside down. In Table 1, "Bridge" refers to the cantilever example.

Chan, Y.-C., Ahmed, F., Wang, L., and Chen, W., 2020, "METASET: An automated data selection method for scalable data-driven design of metamaterials." _ASME 2020 Design Automation Conference_. Paper Number: DETC2020-22681.

    @inproceedings{Chan2020idetc
      year = {2020},
      author = {Yu-Chin Chan and Faez Ahmed and Liwei Wang and Wei Chen},
      title = {METASET: An automated data selection method for scalable data-driven design of metamaterials},
      booktitle={ASME 2020 International Design Engineering Technical Conferences \& Computers and Information in Engineering Conference},
      month={August},
      address={St. Louis, MO, USA}
    }
