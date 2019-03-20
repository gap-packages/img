gap> START_TEST("fr:p1");
gap> SetP1Points(PMCOMPLEX);
gap> SetFloats(PMCOMPLEX);
gap> n := P1Point(0.0);
<0+0i>
gap> s := P1Antipode(n);
<P1infinity>
gap> e := P1Point(1.0);;
gap> w := P1Point(-1.0);;
gap> f := P1Point(1.0i);;
gap> b := P1Point(-1.0i);;
gap> pts := [n,s,e,w,f,b,P1Point(1.e-3),P1Point(1.e-3i),P1Point(1.e-3+1.e-3i),P1Point(1.e-6i)];
[ <0+0i>, <P1infinity>, <1+0i>, <-1+0i>, <0+1i>, <-0-1i>, <0.001+0i>, 
  <0+0.001i>, <0.001+0.001i>, <0+1e-06i> ]
gap> P1Midpoint(s,n);
fail
gap> P1Midpoint(s,e);
<2.4142+0i>
gap> P1Midpoint(n,w);
<-0.41421+0i>
gap> P1Midpoint(pts[9],pts[10]);
<0.0005+0.0005005i>
gap> P1Midpoint(pts[10],pts[10]);
<0+1e-06i>
gap> List(pts,SphereP1);
[ [ 0., 0., 1. ], [ 0., 0., -1. ], [ 1., 0., 0. ], [ -1., 0., 0. ], 
  [ 0., 1., 0. ], [ -0., -1., 0. ], [ 0.002, 0., 0.999998 ], 
  [ 0., 0.002, 0.999998 ], [ 0.002, 0.002, 0.999996 ], [ 0., 2.e-06, 1. ] ]
gap> List(pts,x->P1Distance(x,P1Sphere(SphereP1(x))));
[ 0., 0., 0., 0., 0., 0., 1.32561e-19, 1.32561e-19, 2.42871e-19, 6.67948e-23 ]
gap> DelaunayTriangulation(pts);
<triangulation with 10 vertices, 48 edges and 16 faces>
gap> DelaunayTriangulation(pts,2.0);
<triangulation with 90 vertices, 528 edges and 176 faces>
gap> 
gap> z := Indeterminate(PMCOMPLEX_PSEUDOFIELD,"z");
z
gap> (z^2-1)/(z-1);
z+1._z
gap> STOP_TEST("p1.tst", 10^8);
fr:p1
msecs: 14