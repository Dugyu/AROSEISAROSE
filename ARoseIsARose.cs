using System;
using System.Collections;
using System.Collections.Generic;

using Rhino;
using Rhino.Geometry;

using Grasshopper;
using Grasshopper.Kernel;
using Grasshopper.Kernel.Data;
using Grasshopper.Kernel.Types;



/// <summary>
/// This class will be instantiated on demand by the Script component.
/// </summary>
public class Script_Instance : GH_ScriptInstance
{
#region Utility functions
  /// <summary>Print a String to the [Out] Parameter of the Script component.</summary>
  /// <param name="text">String to print.</param>
  private void Print(string text) { /* Implementation hidden. */ }
  /// <summary>Print a formatted String to the [Out] Parameter of the Script component.</summary>
  /// <param name="format">String format.</param>
  /// <param name="args">Formatting parameters.</param>
  private void Print(string format, params object[] args) { /* Implementation hidden. */ }
  /// <summary>Print useful information about an object instance to the [Out] Parameter of the Script component. </summary>
  /// <param name="obj">Object instance to parse.</param>
  private void Reflect(object obj) { /* Implementation hidden. */ }
  /// <summary>Print the signatures of all the overloads of a specific method to the [Out] Parameter of the Script component. </summary>
  /// <param name="obj">Object instance to parse.</param>
  private void Reflect(object obj, string method_name) { /* Implementation hidden. */ }
#endregion

#region Members
  /// <summary>Gets the current Rhino document.</summary>
  private readonly RhinoDoc RhinoDocument;
  /// <summary>Gets the Grasshopper document that owns this script.</summary>
  private readonly GH_Document GrasshopperDocument;
  /// <summary>Gets the Grasshopper script component that owns this script.</summary>
  private readonly IGH_Component Component;
  /// <summary>
  /// Gets the current iteration count. The first call to RunScript() is associated with Iteration==0.
  /// Any subsequent call within the same solution will increment the Iteration count.
  /// </summary>
  private readonly int Iteration;
#endregion

  /// <summary>
  /// This procedure contains the user code. Input parameters are provided as regular arguments,
  /// Output parameters as ref arguments. You don't have to assign output parameters,
  /// they will have a default value.
  /// </summary>
  private void RunScript(Surface srf, int uCount, int vCount, bool type, int maxIter, ref object CUTLINES, ref object CUTCURVES, ref object TRICENTERS)
  {

    // Final Triangles
    List<SrfTriangle> almostFlatTriangles = new List<SrfTriangle>();

    // Input & Output Triangles
    List<SrfTriangle> inputTriangles = new List<SrfTriangle>();
    List<SrfTriangle> outputTriangles = new List<SrfTriangle>();

    // Initialize Triangles on Surface
    List<SrfTriangle> srfTris = getSrfTriangles(uCount, vCount, getSrfPoints2d(uCount, vCount), type);

    // Iteration
    inputTriangles.AddRange(srfTris);

    for (int i = 0; i < maxIter; ++i)

    {
      //Update Benchmarks
      double[] GauDev = ComputeGaussianDeviation(srf, inputTriangles);
      double gau = GauDev[0];
      double dev = GauDev[1];
      outputTriangles.Clear();
      foreach(SrfTriangle tri in inputTriangles){

        if (tri.Deviation(srf) > dev ){
          SrfTriangle[] triChildren = tri.getChildren();
          outputTriangles.AddRange(triChildren);
        }
        else {
          outputTriangles.Add(tri);
        }
      }

      List<SrfTriangle> swap = inputTriangles;
      inputTriangles = outputTriangles;
      outputTriangles = swap;
    }

    almostFlatTriangles = outputTriangles;
    List <Line> cutLines = new List<Line>();
    List <Curve> cutCurves = new List<Curve>();
    List <Point3d> triCenters = new List<Point3d>();

    foreach (SrfTriangle tri in almostFlatTriangles)
    {
      cutLines.AddRange(tri.SrfTriLines(srf));
      //cutCurves.AddRange(tri.SrfTriCurves(srf));
      triCenters.AddRange(tri.samplePoints(srf));
    }

    CUTLINES = cutLines;
    //CUTCURVES = cutCurves;
    TRICENTERS = triCenters;


  }

  // <Custom additional code> 

  public double[] ComputeGaussianDeviation(Surface srf, List<SrfTriangle> triangles)
  {
    double DevSum = 0.0;
    double GauSum = 0.0;

    foreach(SrfTriangle tri in triangles)
    {
      GauSum += tri.Gaussian(srf);
      DevSum += tri.Deviation(srf);
    }

    double[] result = new double[3];
    result[0] = GauSum * (1.0 / triangles.Count);
    result[1] = DevSum * (1.0 / triangles.Count);

    return result;
  }


  public List<Point2d> getSrfPoints2d(int uCount, int vCount)
  {

    double u0 = 0.0;
    double u1 = 1.0;

    double v0 = 0.0;
    double v1 = 1.0;

    double du = (u1 - u0) / (uCount - 1.0);
    double dv = (v1 - v0) / (vCount - 1.0);

    List<Point2d> points = new List<Point2d>();

    for(int j = 0; j < vCount; ++j) {
      for(int i = 0; i < uCount; ++i) {

        double u = u0 + i * du;
        double v = v0 + j * dv;

        points.Add(new Point2d(u, v));
      }
    }
    return points;
  }


  public List<SrfTriangle> getSrfTriangles(int uCount, int vCount, List<Point2d> pts, bool type)
  {
    List<SrfTriangle> triangles = new List<SrfTriangle>();

    if (type == true){
      for(int j = 0; j < vCount; ++j) {
        for(int i = 0; i < uCount; ++i) {
          if ((i + j) % 2 == 0 && j + 1 < vCount ){
            // Current point
            int A = uCount * j + i;
            int B = A + uCount + 1;
            int C = A + 2;
            int D = A + uCount - 1;
            //Left Corner : 1 Triangle
            if (i - 1 < 0 && i + 2 < uCount){
              triangles.Add(new SrfTriangle(pts[A], pts[B], pts[C]));
            }
              //Majority: 2 Triangles
            else if(i - 1 >= 0 && i + 2 < uCount){
              triangles.Add(new SrfTriangle(pts[A], pts[B], pts[C]));
              triangles.Add(new SrfTriangle(pts[A], pts[D], pts[B]));
            }
              //Right Corner: 1 Triangle
            else if(i - 1 >= 0 && i + 1 < uCount){
              triangles.Add(new SrfTriangle(pts[A], pts[D], pts[B]));
            }
          }
        }
      }
    }

    else{
      for(int j = 0; j < vCount; ++j) {
        for(int i = 0; i < uCount; ++i) {
          if ((i + j) % 2 == 1 && j > 0 ){
            // Current point
            int A = uCount * j + i;
            int B = A - uCount + 1;
            int C = A + 2;
            int D = A - uCount - 1;
            //Left Corner : 1 Triangle
            if (i - 1 < 0 && i + 2 < uCount){
              triangles.Add(new SrfTriangle(pts[A], pts[B], pts[C]));
            }
              //Majority: 2 Triangles
            else if(i - 1 >= 0 && i + 2 < uCount){
              triangles.Add(new SrfTriangle(pts[A], pts[B], pts[C]));
              triangles.Add(new SrfTriangle(pts[A], pts[D], pts[B]));
            }
              //Right Corner: 1 Triangle
            else if(i - 1 >= 0 && i + 1 < uCount){
              triangles.Add(new SrfTriangle(pts[A], pts[D], pts[B]));
            }
          }
        }
      }
    }
    return triangles;
  }


  //Declare triangle on surface class
  public class SrfTriangle {

    // A.B.C Node U V
    // Auto-implemented read&write property:
    public double Au {get;set;}
    public double Av {get;set;}
    public double Bu {get;set;}
    public double Bv {get;set;}
    public double Cu {get;set;}
    public double Cv {get;set;}

    public Point2d A {get;set;}
    public Point2d B {get;set;}
    public Point2d C {get;set;}

    // Constructor that takes no arguments:
    public SrfTriangle()
    {
      Au = 0.0;
      Av = 0.0;
      Bu = 0.0;
      Bv = 0.0;
      Cu = 0.0;
      Cv = 0.0;

    }
    // Constructor that takes 6 arguments:
    public SrfTriangle(double _Au, double _Av, double _Bu, double _Bv, double _Cu, double _Cv)
    {
      Au = _Au;
      Av = _Av;
      Bu = _Bu;
      Bv = _Bv;
      Cu = _Cu;
      Cv = _Cv;
    }
    // Constructor that takes 3 arguments:
    public SrfTriangle(Point2d _A, Point2d _B, Point2d _C)
    {
      Au = _A.X;
      Av = _A.Y;
      Bu = _B.X;
      Bv = _B.Y;
      Cu = _C.X;
      Cv = _C.Y;
    }

    //Extract the points in surface parameter space using Normalized Parameters:
    public Point2d[] SrfTriPoints2d(Surface srf)
    {
      Point2d[] NodePoints2d = new Point2d[3];
      double uMin = srf.Domain(0).Min;
      double uMax = srf.Domain(0).Max;
      double vMin = srf.Domain(1).Min;
      double vMax = srf.Domain(1).Max;

      double nAu = (uMax - uMin) * Au + uMin;
      double nAv = (vMax - vMin) * Av + vMin;
      double nBu = (uMax - uMin) * Bu + uMin;
      double nBv = (vMax - vMin) * Bv + vMin;
      double nCu = (uMax - uMin) * Cu + uMin;
      double nCv = (vMax - vMin) * Cv + vMin;

      NodePoints2d[0] = new Point2d(nAu, nAv);
      NodePoints2d[1] = new Point2d(nBu, nBv);
      NodePoints2d[2] = new Point2d(nCu, nCv);

      return NodePoints2d;
    }

    //Extract the actual Point3d of SrfTriangle Nodes using Normalized Parameters:
    public Point3d[] SrfTriPoints(Surface srf)
    {
      Point3d[] Nodes = new Point3d[3];
      Point2d[] Nodes2d = this.SrfTriPoints2d(srf);

      Point2d A = Nodes2d[0];
      Point2d B = Nodes2d[1];
      Point2d C = Nodes2d[2];

      Nodes[0] = srf.PointAt(A.X, A.Y);
      Nodes[1] = srf.PointAt(B.X, B.Y);
      Nodes[2] = srf.PointAt(C.X, C.Y);

      return Nodes;
    }

    //Compute the Curve Edges of SrfTriangle:
    public Curve[] SrfTriCurves(Surface srf)
    {
      Curve[] SrfTriCurves = new Curve[3];
      Point2d[] Nodes2d = this.SrfTriPoints2d(srf);

      Point2d A = Nodes2d[0];
      Point2d B = Nodes2d[1];
      Point2d C = Nodes2d[2];

      SrfTriCurves[0] = srf.ShortPath(A, B, 0.001);
      SrfTriCurves[1] = srf.ShortPath(B, C, 0.001);
      SrfTriCurves[2] = srf.ShortPath(C, A, 0.001);

      return SrfTriCurves;
    }


    //Compute the straight Line Edges of SrfTriangle:
    public Line[] SrfTriLines(Surface srf)
    {

      Line[] SrfTriLines = new Line[3];
      Point3d[] Nodes = this.SrfTriPoints(srf);

      SrfTriLines[0] = new Line(Nodes[0], Nodes[1]);
      SrfTriLines[1] = new Line(Nodes[1], Nodes[2]);
      SrfTriLines[2] = new Line(Nodes[2], Nodes[0]);

      return SrfTriLines;
    }

    //Compute Central Point2d of a triangular on a srf
    public Point2d SrfTriCenter2d(Surface srf)

    {
      Point2d[] Nodes2d = this.SrfTriPoints2d(srf);
      Point2d center2d = (Nodes2d[0] + Nodes2d[1] + Nodes2d[2]) * (1 / 3.0);
      return center2d;
    }

    //Compute Central Point3d of a triangular on a srf
    public Point3d SrfTriCenter(Surface srf)
    {
      Point2d center2d = this.SrfTriCenter2d(srf);
      Point3d center = srf.PointAt(center2d.X, center2d.Y);
      return center;
    }

    public Point2d[] samplePoints2d(Surface srf)
    {
      Point2d[] Nodes = this.SrfTriPoints2d(srf);

      Point2d midAB = (Nodes[0] + Nodes[1]) * 0.5;
      Point2d midBC = (Nodes[1] + Nodes[2]) * 0.5;
      Point2d midCA = (Nodes[2] + Nodes[0]) * 0.5;

      Point2d nA = (Nodes[0] + midAB + midCA ) * (1.0 / 3.0);
      Point2d nB = (Nodes[1] + midAB + midBC ) * (1.0 / 3.0);
      Point2d nC = (Nodes[2] + midCA + midBC ) * (1.0 / 3.0);

      Point2d[] sample2d = new Point2d[7];
      sample2d[0] = this.SrfTriCenter2d(srf);
      sample2d[1] = midAB;
      sample2d[2] = midBC;
      sample2d[3] = midCA;
      sample2d[4] = nA;
      sample2d[5] = nB;
      sample2d[6] = nC;
      return sample2d;
    }

    public Point3d[] samplePoints(Surface srf)
    {
      Point2d[] sample2d = this.samplePoints2d(srf);
      Point3d[] sample = new Point3d[sample2d.Length];
      for(int i = 0, max = sample2d.Length; i < max; i++) {
        sample[i] = srf.PointAt(sample2d[i].X, sample2d[i].Y);
      }
      return sample;
    }

    //Compute Relative Perimeter Length

    public double perimeter2d()
    {
      return (Au - Bu) * (Au - Bu) + (Av - Bv) * (Av - Bv) +
        (Bu - Cu) * (Bu - Cu) + (Bv - Cv) * (Bv - Cv) +
        (Au - Cu) * (Au - Cu) + (Av - Cv) * (Av - Cv);
    }


    public double perimeter(Surface srf)
    {
      Point3d[] Nodes = this.SrfTriPoints(srf);
      double dist = (Nodes[0].X - Nodes[1].X) * (Nodes[0].X - Nodes[1].X) + (Nodes[0].Y - Nodes[1].Y) * (Nodes[0].Y - Nodes[1].Y)
        + (Nodes[1].X - Nodes[2].X) * (Nodes[1].X - Nodes[2].X) + (Nodes[1].Y - Nodes[2].Y) * (Nodes[1].Y - Nodes[2].Y)
        + (Nodes[0].X - Nodes[2].X) * (Nodes[0].X - Nodes[2].X) + (Nodes[0].Y - Nodes[2].Y) * (Nodes[0].Y - Nodes[2].Y);
      return dist;

    }


    //Compute Relative Gaussian Curvature at the Central Point
    public double Gaussian(Surface srf)
    {
      Point2d[] sample2d = this.samplePoints2d(srf);
      double gauSum = 0.0;
      foreach (Point2d pt in sample2d)
      {
        SurfaceCurvature cv = srf.CurvatureAt(pt.X, pt.Y);
        gauSum += Math.Abs(cv.Gaussian);
      }
      gauSum *= (1.0 / this.perimeter(srf));
      return gauSum;
    }

    //Compute the Deviation (Dist to its triangular plane) at the Central Point
    public double Deviation(Surface srf)
    {

      Point3d[] sample = this.samplePoints(srf);
      Point3d[] Nodes = this.SrfTriPoints(srf);
      Plane triPlane = new Plane(Nodes[0], Nodes[1], Nodes[2]);

      double devSum = 0.0;

      foreach (Point3d pt in sample)
      {
        devSum += Math.Abs(triPlane.DistanceTo(this.SrfTriCenter(srf)));
      }
      devSum *= (1.0 / this.perimeter(srf));
      return devSum;
    }

    //Produce Child SrfTriangles
    public SrfTriangle[] getChildren()
    {
      double midABu = (Au + Bu) * 0.5;
      double midABv = (Av + Bv) * 0.5;
      double midBCu = (Bu + Cu) * 0.5;
      double midBCv = (Bv + Cv) * 0.5;
      double midCAu = (Cu + Au) * 0.5;
      double midCAv = (Cv + Av) * 0.5;

      SrfTriangle[] childs = new SrfTriangle[4];
      childs[0] = new SrfTriangle(midABu, midABv, midBCu, midBCv, midCAu, midCAv);
      childs[1] = new SrfTriangle(Au, Av, midABu, midABv, midCAu, midCAv);
      childs[2] = new SrfTriangle(midABu, midABv, Bu, Bv, midBCu, midBCv);
      childs[3] = new SrfTriangle(midCAu, midCAv, midBCu, midBCv, Cu, Cv);

      return childs;
    }

  }

  // </Custom additional code> 
}