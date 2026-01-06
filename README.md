A Python implementation of the algorithm developed by Oliver Jenkinson and Mark Pollicott in their 2001 paper<sup>\[1\]</sup>, written as part of my MSc Pure Mathematics dissertation<sup>\[2\]</sup>, supervised by Charles Walkden.

For arbitrary-precision floating-point arithmetic, [mpmath](https://github.com/mpmath/mpmath) is used. The code is thoroughly commented, and a more detailed breakdown of the implementation is given in section 4.2 of my dissertation<sup>\[2\]</sup>. For the set E<sub>2</sub>, this implementation agrees with the results in both the original paper by Jenkinson & Pollicott<sup>\[1\]</sup>, and their more recent paper from 2018<sup>\[3\]</sup>.

\[1\] Jenkinson, O. & Pollicott, M. (2001). Computing the dimension of dynamically defined sets: E<sub>2</sub> and bounded continued fractions. Ergodic Theory and Dynamical Systems, 21(5), 1429–1445. doi:10.1017/S0143385701001687

\[2\] Woodell, R. (2023). A computational approach to finding the Hausdorff dimension of continued fraction Cantor sets, ([available in this repository](https://github.com/rowanwoodell/jenkinson-pollicott/blob/main/MSc_Dissertation_FINAL.pdf))

\[3\] Jenkinson, O. & Pollicott, M. (2017). Rigorous effective bounds on the Hausdorff dimension of continued fraction Cantor sets: a hundred decimal digits for the dimension of E<sub>2</sub>, [arXiv:1611.09276](https://arxiv.org/abs/1611.09276)

