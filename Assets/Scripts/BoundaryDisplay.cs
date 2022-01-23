using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class BoundaryDisplay : MonoBehaviour
{
    public Boundary boundary;
    public Color    color = Color.green;
    public int      vertexCount;

    public void Clear()
    {
        boundary = null;
    }

    public void OnDrawGizmos()
    {
        if (boundary == null) return;

        Gizmos.color = color;

        var prevMatrix = Gizmos.matrix;
        Gizmos.matrix = transform.localToWorldMatrix;

        vertexCount = 0;

        for (int i = 0; i < boundary.Count; i++)
        {
            var polyline = boundary.Get(i);
            int count = polyline.Count;
            for (int j = 0; j < count; j++)
            {
                Gizmos.DrawLine(polyline[j], polyline[(j + 1) % count]);

                vertexCount++;
            }
        }

        Gizmos.matrix = prevMatrix;
    }
}
