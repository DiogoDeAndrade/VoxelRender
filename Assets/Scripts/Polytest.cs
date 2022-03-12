using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using NaughtyAttributes;

public class Polytest : MonoBehaviour
{
    public bool isCW;

    [Button("Compute winding")]
    void ComputeWinding()
    {
        PolylineDisplay pd = GetComponent<PolylineDisplay>();
        var polyline = pd.polyline;

        isCW = polyline.isCW();
    }

    [Button("Reverse winding")]
    void ReverseWinding()
    {
        PolylineDisplay pd = GetComponent<PolylineDisplay>();
        var polyline = pd.polyline;

        polyline.ReverseOrder();

        isCW = polyline.isCW();
    }

    [Button("Build Mesh")]
    void BuildMesh()
    {
        PolylineDisplay pd = GetComponent<PolylineDisplay>();
        var polyline = pd.polyline;
        var mesh = MeshTools.TriangulateEarClipping(polyline);

        GetComponent<MeshFilter>().sharedMesh = mesh;
    }
}
