# Ende: Motion Planning for a Magnetorquer-actuated Satellite

Application of kynodynamical motion planning algorithms to perform attitude maneuvers for a magnetorquer-actuated LEO satellite with pointing constrains. [The full paper is available here.](paper.pdf)


## Description
<img src="ende_tree.png" alt="A tree illustrated by Ende." width="200" align="right"/>

Demonstrates a novel path planning algorithm, Bonsai-RRT, to solving the magnetorquer-actuated LEO satellite problem.

For single runs, execute `experiments.jl` with `alg_mode` and `control` mode set accordingly. Other parameters can also be changed here.

For benchmark runs, execute `benchmark.jl` with cases defined in the `cases` dictionary. `benchmark.jl` includes an example of the `bonsai-pd` case setup.

Named after the 10c. illuminator/cartographer Ende.


## Authors

Developed by Mark Stephenson as a final project for CU Boulder's ASEN 5254: Algorithmic Motion Planning in Fall of 2022, taught by Morteza Lahijanian.

## References
Incomplete list, see paper for more.

### Constained Attitude Maneuvers
* R. Calaon and H. Schaub, “Constrained Attitude Maneuvering via Modified-Rodrigues-Parameter-Based Motion Planning Algorithms,” Journal of Spacecraft and Rockets, vol. 59, no. 4, pp. 1342–1356, Jul. 2022, doi: 10.2514/1.A35294.

### Mangentorquer Control and Dynamics
* F. L. Markley and J. L. Crassidis, Fundamentals of Spacecraft Attitude Determination and Control. New York, NY: Springer New York, 2014. doi: 10.1007/978-1-4939-0802-8.
* L. Musser and W. L. Ebert, “Autonomous Spacecraft Attitude Control Using Magnetic Torquing Only,” p. 16, 1989.
* A. C. Stickler and K. T. Alfriend, “Elementary Magnetic Attitude Control System,” Journal of Spacecraft and Rockets, vol. 13, no. 5, pp. 282–287, May 1976, doi: 10.2514/3.57089.
* M. F. Erturk and C. Hajiyev, “Magnetic Torquers Only Attitude Control of a 3U Cube Satellite,” WSEAS TRANSACTIONS ON SIGNAL PROCESSING, vol. 18, pp. 128–133, Jun. 2022, doi: 10.37394/232014.2022.18.18.

### Kinodynamic Planning
* Y. Li, Z. Littlefield, and K. E. Bekris, “Asymptotically Optimal Sampling-based Kinodynamic Planning.” arXiv, Feb. 05, 2016. Accessed: Oct. 31, 2022. [Online]. Available: http://arxiv.org/abs/1407.2896
