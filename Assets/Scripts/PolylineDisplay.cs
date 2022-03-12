using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class PolylineDisplay : MonoBehaviour
{
    public Polyline polyline;
    public Color    color = Color.green;

    public void Clear()
    {
        polyline = null;
    }

    public void OnDrawGizmos()
    {
        if (polyline == null) return;

        Gizmos.color = color;

        var prevMatrix = Gizmos.matrix;
        Gizmos.matrix = transform.localToWorldMatrix;

        int count = polyline.Count;
        for (int j = 0; j < count; j++)
        {
            Gizmos.DrawLine(polyline[j], polyline[(j + 1) % count]);
        }

        Gizmos.matrix = prevMatrix;
    }
}
