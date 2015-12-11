To use this you will need sage and the latest version of mine and Dylan's  ClusterAlgebra package. Get it while it is hot at 

https://raw.githubusercontent.com/Etn40ff/cluster_seed_reborn/cluster_algebras/cluster_algebra.py

and remember to change lines 413-416.

To use these routines:

$ sage
sage: %attach /path/to/cluster_algebra.py
sage: %attach /path/to/find_sortable_cones.py
sage: B = Matrix([[0,-1,1],[1,0,-1],[-2,2,0]])
sage: S = SortableCones(B)

profit.

BIG WARNING: this is experimental code ant it is SLOWWWWWWW!

