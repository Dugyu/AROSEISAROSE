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
  private void RunScript(Surface srf, List<Point3d> pts, int n, double roseRange, List<Vector3d> normals, ref object CONNECTIONS, ref object ROSENODES, ref object ROSETFORM, ref object ROSESCALE)
  {

    List<List<int>> neighbors = new List<List<int>> ();
    List<List<double>> dists = new List<List<double>> ();

    // Traverse point list to select point center
    for (int j = 0; j < pts.Count; j++){

      // Store the index, distance of each neighbor to current center point
      List <int> neighborsLocal = new List<int>();
      List <double> distsLocal = new List<double>();

      for (int i = 0; i < pts.Count; i++){
        // not calculating itself
        if(i != j)
        {
          double dist = pts[j].DistanceTo(pts[i]);
          neighborUpdate(i, dist, n, ref neighborsLocal, ref distsLocal);
        }
      }

      // Add the Local Neighbors to the Global Neighbors List
      neighbors.Add(neighborsLocal);
      dists.Add(distsLocal);
    }
    // Clean Neighbors, Ensure no more than N lines is connected to each point
    //neighborClean(n, ref neighbors);

    // output Rose Nodes
    List <Point3d> roseNodes = new List<Point3d>();
    // output neighbor lines
    List<Line> lines = new List<Line>();
    // output Roses
    List<Transform> roseTransforms = new List<Transform>();
    List<double> roseScale = new List<double>();

    for (int j = 0; j < pts.Count; j++){
      Line[] neighborLines = new Line[n];
      double sum = 0;
      for (int i = 0; i < neighbors[j].Count; i++){
        Line link = new Line(pts[j], pts[neighbors[j][i]]);
        sum += link.Length;
        neighborLines[i] = link;
      }
      lines.AddRange(neighborLines);

      if (sum > roseRange ){
        roseNodes.Add(pts[j]);

        Plane frame;
        frame = new Plane(pts[j], normals[j]);
        //srf.FrameAt(ptsUV[j].X, ptsUV[j].Y, out frame);

        Transform tform = Transform.PlaneToPlane(Plane.WorldXY, frame);
        roseTransforms.Add(tform);
        roseScale.Add(sum * 0.005);
      }
    }

    CONNECTIONS = lines;
    ROSENODES = roseNodes;
    ROSETFORM = roseTransforms;
    ROSESCALE = roseScale;

  }

  // <Custom additional code> 
  public void neighborClean(int n, ref List<List<int>> neighbors)
  {
    //Ensure Every Point is only connected to equal to N neighbors or less

    for (int j = 0; j < neighbors.Count; j++){
      for (int i = n - 1; i >= 0; i--){

        // Weather The neighbor current point links to Contain the current point or Not
        bool isContained = false;
        for (int k = 0; k < neighbors[neighbors[j][i]].Count; k++)
        {
          if (neighbors[neighbors[j][i]][k] == j){
            isContained = true;
            break;
          }
        }

        // the neighbor that the current point links to does not contain current point
        if (!isContained)
        {
          // Not Full
          if(neighbors[neighbors[j][i]].Count < n)
          {
            neighbors[neighbors[j][i]].Add(j);
          }
            // Full
          else
          {
            neighbors[j].RemoveAt(i);
          }
        }

      }
    }

  }


  public void neighborUpdate(int index, double dist, int n, ref List<int> neighbors, ref List<double> dists)
  {
    // Empty Neighbors
    if (dists.Count == 0){
      dists.Add(dist);
      neighbors.Add(index);
      return;
    }
      // Part Empty Neighbors
    else if(dists.Count < n){
      //If Closer than the farest neighbor, starting the sort, find its place in the neighbors list
      if(dist < dists[dists.Count - 1])
      {
        for (int i = 0; i < dists.Count; i++)
        {
          // Insert according the distance order
          if (dist < dists[i])
          {
            dists.Insert(i, dist);
            neighbors.Insert(i, index);
            break;
          }
        }
      }
        //Add to the Farest Point
      else{
        dists.Add(dist);
        neighbors.Add(index);
      }
    }


      // Full Neighbors
      //If Closer than the farest neighbor, starting the sort, find its place in the neighbors list
    else if(dist < dists[dists.Count - 1])
    {

      // Traverse neighbors
      for (int i = 0; i < dists.Count; i++)
      {
        // Insert according the distance order
        if (dist < dists[i])
        {
          dists.Insert(i, dist);
          neighbors.Insert(i, index);
          break;
        }
      }
      // Remove the farest point
      dists.RemoveAt(dists.Count - 1);
      neighbors.RemoveAt(neighbors.Count - 1);
    }
  }

  // </Custom additional code> 
}