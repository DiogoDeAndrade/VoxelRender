using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class LineIntersectionTest : MonoBehaviour
{
    public Vector3 p1;
    public Vector3 p2;
    public Vector3 p3;
    public Vector3 p4;

    private void OnDrawGizmosSelected()
    {
        Vector3 intersection;
        if (Line.Intersect(p1, p2, p3, p4, out intersection))
        {
            Gizmos.color = Color.yellow;
            Gizmos.DrawSphere(intersection, 0.01f);

            Gizmos.color = Color.red;
        }
        else
        {
            Gizmos.color = Color.green;
        }

        Gizmos.DrawLine(p1, p2);
        Gizmos.DrawLine(p3, p4);
    }
}
