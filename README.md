# Computation of the Hausdorff Distance between Two Compact Convex Sets

The Hausdorff distance between two closed sets has important theoretical and practical applications.
Yet apart from finite point clouds, there appear to be no generic algorithms for computing this quantity.
Because many infinite sets are defined by algebraic equalities and inequalities, this is a huge gap.
The current paper constructs Frank-Wolfe and projected gradient ascent algorithms for computing the Hausdorff distance between two compact convex sets.
Although these algorithms are guaranteed to go uphill, they can get trapped by local maxima.
To avoid this defect, we investigate a homotopy method that gradually deforms two balls into the two target sets.
The Frank-Wolfe and projected gradient algorithms are tested on two pairs $(A, B)$ of compact convex sets, where: (1) $A$ is the box $[−1, 1]$ translated by 1, and $B$ is the intersection of the unit ball and the nonnegative orthant; and (2) $A$ is the probability simplex and $B$ is the $\ell_{1}$ unit ball translated by 1.
For problem (2) we find the Hausdorff distance analytically. Projected gradient ascent is more reliable than Frank-Wolfe and finds the exact solution of problem (2).
Homotopy improves the performance of both algorithms when the exact solution is unknown or unattained.
