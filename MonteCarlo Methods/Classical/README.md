+ #### 2DIsingModelPBC.f90
simplest example of how monte carlo method is applied on 2D Ising model with periodic boundary condition using metropolis sampling. The peak of specific heat $c$ and magnetic susceptibility $\chi$, you can find the critical point $T_c$. The intersection point of [Binder cumulant](http://www.sklogwiki.org/SklogWiki/index.php/Binder_cumulant) curves under different system sizes $L$ will give exact critical point $T_c$.

+ #### 2DIsingSW.f90
an implementation of [Swendsen-Wang Algorithm](http://en.wikipedia.org/wiki/Swendsen%E2%80%93Wang_algorithm), [paper](http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.58.86)

+ #### 2DHeisenbergWolffPBC.f90
an implementation of [Wolff Algorithm](http://en.wikipedia.org/wiki/Wolff_algorithm), [paper](http://journals.aps.org/prl/abstract/10.1103/PhysRevLett.62.361)

+ #### WolffMixMetropolis.f90
for an anisotropic 2D Heisenberg model, this code shows a hybrid method using both metropolis local update ($z$ component) and wolff cluster update ($x-y$ component).
