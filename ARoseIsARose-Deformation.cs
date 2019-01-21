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
  private void RunScript(Surface srf, int uCount, int vCount, double torque, double offset, double focus, double force, ref object OFFSETpts, ref object OFFSETsrf, ref object UCOUNT, ref object VCOUNT, ref object PTS2D)
  {
    double u0 = srf.Domain(0).Min;
    double u1 = srf.Domain(0).Max;

    double v0 = srf.Domain(1).Min;
    double v1 = srf.Domain(1).Max;

    double du = (u1 - u0) / (uCount - 1.0);
    double dv = (v1 - v0) / (vCount - 1.0);

    List<Point3d> points = new List<Point3d>();
    List<Vector3d> normals = new List<Vector3d>();
    List<Point3d> offsetPoints = new List<Point3d>();
    List<Point2d> points2d = new List<Point2d>();

    for (int j = 0; j < vCount; ++j)
    {
      for (int i = 0; i < uCount; ++i)
      {

        double u = u0 + i * du;
        double v = v0 + j * dv;

        // transfer to interval (0,1)
        double uu = srf.Domain(0).NormalizedParameterAt(u);
        double vv = srf.Domain(1).NormalizedParameterAt(v);

        // find the focus point along midline according to user input
        double nu = uu - 0.5;
        double nv = vv - focus;

        Point3d p = srf.PointAt(u, v);
        Point2d pp = new Point2d(uu, vv);

        Vector3d n = srf.NormalAt(u, v);
        SurfaceCurvature cv = srf.CurvatureAt(u, v);

        // find the principal curvature directions at current uv
        Vector3d k1dir = cv.Direction(0);
        Vector3d k2dir = cv.Direction(1);

        // apply torque to points along principal curvature directions
        /* make it symmetric
        if uCount is odd, make the midline points static*/
        if (uCount % 2 == 1)
        {
          if (i != uCount / 2)
          {
            p += (k2dir + k1dir) * torque;
          }
        }
          // if uCount is even, make the two lines of points in the middle static
        else
        {
          if (i != uCount / 2 && i != (uCount / 2) - 1)
          {
            p += (k2dir + k1dir) * torque;
          }
        }

        // apply force to offset with regard to the focus point
        double dd = nu * nu + nv * nv;
        double off = (0.1 + Math.Exp(-force * dd)) * offset;

        points2d.Add(pp);
        points.Add(p);
        offsetPoints.Add(p + n * (-cv.Mean) * off); //apply mean curvature flow with offset
      }
    }


    Surface offsetSurface = NurbsSurface.CreateThroughPoints(offsetPoints, vCount, uCount, 3, 3, false, false);

    OFFSETpts = offsetPoints;
    OFFSETsrf = offsetSurface;
    UCOUNT = uCount;
    VCOUNT = vCount;
    PTS2D = points2d;

  }

  // <Custom additional code> 

  // </Custom additional code> 
}