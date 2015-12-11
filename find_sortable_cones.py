class SortableCones(SageObject):
    def __init__(self, B):
        self._B = B
        self._D = diagonal_matrix(B.is_skew_symmetrizable(return_diag=True))
        self._mBT = -B.transpose()
        self._DT = diagonal_matrix(self._mBT.is_skew_symmetrizable(return_diag=True))
        A = ClusterAlgebra(B)
        self._g_vector_fan = A.cluster_fan()
        A_dual = ClusterAlgebra(-B.transpose())
        c_vectors_dual = []
        for S in A_dual.seeds():
            c_vectors_dual += S.c_vectors()
        c_vectors_dual = set([ c for c in c_vectors_dual if all([cc>=0 for cc in c]) ])
        self._c_vectors_dual = c_vectors_dual
        self._n = B.ncols()
        hyperplene_arrangement_space = HyperplaneArrangements(QQ, var(['x%d'%i for i in range(self._n)]) )
        hyperplane_arrangement = hyperplene_arrangement_space([(c,0) for c in self._c_vectors_dual])
        self._supporting_fan = Fan(map(lambda C: Cone(C),  [ R.rays() for R in hyperplane_arrangement.regions()] ))
        self._oriented_triples = []
        var('x,y')
        for (a,b,c) in tuples(self._c_vectors_dual,3):
            if len(set([a,b,c])) == 3:
                if vector(a)*self._DT*self._mBT*vector(c) > 0 and solve(map(lambda rel: rel==0, vector(a)*x+vector(c)*y-vector(b))+[x>0,y>0],x,y) != []:
                    self._oriented_triples.append((a,b,c))
              
    def number_of_sortable(self):
        return len(self.sortable_cones())

    @cached_method
    def sortable_cones(self):
        return tuple([ C for C in self._supporting_fan.cones()[-1] if self._is_sortable(C) ])
    
    @cached_method
    def killing_triples(self, cone):
        inversions = self._inversions(cone)
        killing_triples = []
        for (a,b,c) in self._oriented_triples:
            if b in inversions and a not in inversions:
                killing_triples.append((a,b,c))
        return tuple(killing_triples)
        
    @cached_method
    def _inversions(self, C):
        return tuple([c for c in self._c_vectors_dual if all([vector(c)*vector(r)<=0 for r in C.rays()])])

    @cached_method
    def _is_sortable(self,cone):
        inversions = self._inversions(cone)
        for (a,b,c) in self._oriented_triples :
            if b in inversions and a not in inversions:
                return False
        return True

    @cached_method
    def is_minimal(self, cone, cones_tuple):
        inversions = set(self._inversions(cone))
        return all([inversions.issubset(set(self._inversions(other_cone))) for other_cone in cones_tuple])
    
    @cached_method
    def projection_fibers(self):
        return tuple([[ C for C in self._supporting_fan.cones()[-1] if self._g_vector_fan.cone_containing(C) == D] for D in self._g_vector_fan.cones()[-1] ])
    
    @cached_method
    def pi_down_image(self):
        return tuple(map(self.minimal_element, self.projection_fibers()))

    @cached_method
    def minimal_element(self, cones_tuple):
        minimal = []
        for cone in cones_tuple:
            if self.is_minimal(cone, cones_tuple):
                minimal.append(cone)
        if len(minimal)!= 1:
            raise ValueError("This does not have a minimum!")
        return minimal[0]
    
    @cached_method
    def minimal_element_faster(self, cones_tuple):
        inversions_tuple = map(lambda Y: set(self._inversions(Y)), cones_tuple)
        min_inversions = inversions_tuple[0]
        for X in inversions_tuple[1:]:
            min_inversions = min_inversions.intersection(X)
        for C in cones_tuple:
            if set(self._inversions(C)) == min_inversions:
                return C
        print cones_tuple
        raise ValueError("This does not have a minimum")


     




################
# Stereographic projections routines
################

def _stereo_coordinates(x, north=(1,0,0), right=(0,1,0), translation=-1):
    r"""
    Project stereographically points from a sphere
    """
    from sage.misc.functional import norm
    north=_normalize(north)
    right=vector(right)
    right=_normalize(right-(right*north)*north)
    if norm(right) == 0:
        raise ValueError ("Right must not be linearly dependent from north")
    top=north.cross_product(right)
    x=_normalize(x)
    p=(translation-north*x)/(1-north*x)*(north-x)+x
    return vector((right*p, top*p ))
 
def _normalize(x):
    r"""
    make x of length 1
    """
    from sage.misc.functional import norm
    x=vector(RR, x)
    if norm(x) == 0:
        return x
    return vector(x/norm(x))

def rays(fan):
    return map(lambda X: X.rays()[0], fan.cones()[1])

def walls(fan):
    return map(lambda X: X.rays(), fan.cones()[2]) 

def plot_fan_stereographically(rays, walls, northsign=1, north=vector((-1,-1,-1)), right=vector((1,0,0)), colors=None, thickness=None):
    from sage.plot.graphics import Graphics
    from sage.plot.point import point
    from sage.misc.flatten import flatten
    from sage.plot.line import line
    from sage.misc.functional import norm
    
    if colors == None:
        colors = dict([('walls','black'),('rays','red')])

    if thickness == None:
        thickness = dict([('walls',0.5),('rays',20)])


    G = Graphics()
    
    for (u,v) in walls:
        G += _stereo_arc(vector(u),vector(v),vector(u+v),north=northsign*north,right=right,color=colors['walls'],thickness=thickness['walls'],zorder=len(G))
   
    for v in rays: 
        G += point(_stereo_coordinates(vector(v),north=northsign*north,right=right),color=colors['rays'],zorder=len(G),size=thickness['rays'])
    
    G.set_aspect_ratio(1)
    G._show_axes = False
    return G



def _arc3d((first_point,second_point),center=(0,0,0),**kwds):
    # FIXME: refactor this before publishing
    r"""
    Usese parametric_plot3d() to plot arcs of a circle.
    We only plot arcs spanning algles strictly less than pi.
    """
    # For sanity purposes convert the input to vectors
    from sage.misc.functional import norm
    from sage.modules.free_module_element import vector
    from sage.functions.trig import arccos, cos, sin
    from sage.plot.plot3d.parametric_plot3d import parametric_plot3d
    center=vector(center)
    first_point=vector(first_point)
    second_point=vector(second_point)
    first_vector=first_point-center
    second_vector=second_point-center
    radius=norm(first_vector)
    if norm(second_vector)!=radius:
        raise ValueError("Ellipse not implemented")
    first_unit_vector=first_vector/radius
    second_unit_vector=second_vector/radius
    normal_vector=second_vector-(second_vector*first_unit_vector)*first_unit_vector
    if norm(normal_vector)==0:
        print (first_point,second_point)
        return
    normal_unit_vector=normal_vector/norm(normal_vector)
    scalar_product=first_unit_vector*second_unit_vector
    if abs(scalar_product) == 1:
        raise ValueError("The points are alligned")
    angle=arccos(scalar_product)
    var('t')
    return parametric_plot3d(center+first_vector*cos(t)+radius*normal_unit_vector*sin(t),(0,angle),**kwds)
     

     
def _arc(p,q,s,**kwds):
    #rewrite this to use polar_plot and get points to do filled triangles
    from sage.misc.functional import det
    from sage.plot.line import line
    from sage.misc.functional import norm
    from sage.symbolic.all import pi
    from sage.plot.arc import arc
     
    p,q,s = map( lambda x: vector(x), [p,q,s])
     
    # to avoid running into division by 0 we set to be colinear vectors that are
    # almost colinear
    if abs(det(matrix([p-s,q-s])))<0.01:
        return line((p,q),**kwds)
     
    (cx,cy)=var('cx','cy')
    equations=[
            2*cx*(s[0]-p[0])+2*cy*(s[1]-p[1]) == s[0]**2+s[1]**2-p[0]**2-p[1]**2,
            2*cx*(s[0]-q[0])+2*cy*(s[1]-q[1]) == s[0]**2+s[1]**2-q[0]**2-q[1]**2
            ]
    c = vector( [solve( equations, (cx,cy), solution_dict=True )[0][i] for i in [cx,cy]] )
     
    r = norm(p-c)
     
    a_p,a_q,a_s = map( _to_angle, [p-c,q-c,s-c])
    angles = [a_p,a_q,a_s]
    angles.sort()
     
    if a_s == angles[0]:
        return arc( c, r, angle=angles[2], sector=(0,2*pi-angles[2]+angles[1]), **kwds)
    if a_s == angles[1]:
        return arc( c, r, angle=angles[0], sector=(0,angles[2]-angles[0]), **kwds)
    if a_s == angles[2]:
        return arc( c, r, angle=angles[1], sector=(0,2*pi-angles[1]+angles[0]), **kwds)
     
def _to_angle((x,y)):
    from sage.functions.trig import arctan, arccot
    from sage.symbolic.all import pi
    if x >= -y and x >= y:
        return arctan(y/x)
    if x >= -y and x < y:
        return arccot(x/y)
    if x < -y and x < y:
        return pi+arctan(y/x)
    if x < -y and x >= y:
        return pi+arccot(x/y)
     
     
     
def _stereo_arc(x,y, xy=None,  north=(1,0,0), right=(0,1,0), translation=-1, **kwds):
    from sage.misc.functional import n
    x=vector(x)
    y=vector(y)
    sx=n(_stereo_coordinates(x, north=north, right=right, translation=translation))
    sy=n(_stereo_coordinates(y, north=north, right=right, translation=translation))
    if xy == None:
        xy=x+y
    sxy=n(_stereo_coordinates(xy, north=north, right=right, translation=translation))
    return _arc(sx,sy,sxy,**kwds)                                                               
