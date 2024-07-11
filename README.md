# MATPOWER-SICNM

 This repo is an extension package for the algorithms of Matpower, includeing:

 1. Our proposed semi-implicit continuous Newton method (SICNM).
 2. The classic Iwamoto's robust method [^iwamoto];
 3. The original explicit[^cnm1] and implicit[^cnm2] Continuous Newton method proposed by Federico Milano .
 4. Other algorithms used for comparison in our paper.

 This repo also holds the source codes of our paper currently under peer review.

## Usage 

Please add the `src` folder to your matlab path. Use the `runpfa()` function to run the cases in `cases` folder.


## Cite Us
[1] R. Yu, W. Gu, S. Lu, and Y. Xu, “Semi-implicit continuous newton method for power flow analysis,” 2023, arXiv:2312.02809.


[^iwamoto]: S. Iwamoto and Y. Tamura, “A load flow calculation method for ill- conditioned power systems,” IEEE Transactions on Power Apparatus and Systems, vol. PAS-100, no. 4, pp. 1736–1743, 1981.
[^cnm1]: F. Milano, “Continuous newton’s method for power flow analysis,” IEEE Trans. Power Syst., vol. 24, no. 1, pp. 50–57, 2009.
[^cnm2]: F. Milano, “Implicit continuous newton method for power flow analysis,” IEEE Trans. Power Syst., vol. 34, no. 4, pp. 3309–3311, 2019.
