### Dynamic Service Migration in Mobile Edge Computing Based on Markov Decision Process

This is the simulation code for the paper S. Wang, R. Urgaonkar, M. Zafer, T. He, K. Chan, K. K. Leung, "Dynamic service migration in mobile edge computing based on Markov decision process," IEEE/ACM Transactions on Networking, vol. 27, no. 3, pp. 1272 – 1288, Jun. 2019. (arXiv link: [https://arxiv.org/abs/1506.05261](https://arxiv.org/abs/1506.05261))

The code runs best on MATLAB. It also runs on [GNU Octave](https://www.gnu.org/software/octave/) but the plot of instantaneous cost may be shown with lower granularity.

To reproduce the random walk result (Fig. 6 of the paper), run ```mainRandomWalk.m```.

To reproduce the result with real base station locations (Fig. 8 of the paper), run ```mainRealCellLocation.m```.

There are certain parameters in ```mainRandomWalk.m``` and ```mainRealCellLocation.m``` that can be changed for different experiments. The main algorithms are implemented in ```algorithms.m``` which is called within ```mainRandomWalk.m``` and ```mainRealCellLocation.m```.

Real user traces are obtained from [http://crawdad.org/epfl/mobility/20090224/](http://crawdad.org/epfl/mobility/20090224/) and the base station locations are obtained from [http://www.antennasearch.com/](http://www.antennasearch.com/). These are saved in ```traceRealCellLocations.mat```.
