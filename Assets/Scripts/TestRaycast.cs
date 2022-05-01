using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using NaughtyAttributes;

public class TestRaycast : MonoBehaviour
{
    public MeshFilter   meshFilter;
    public float        maxDistance = 100.0f;

    public int submeshId;
    public int triId;

    MeshOctree  meshOctree;

    [Button("Rebuild Octree")]
    void RebuildOctree()
    {
        meshOctree = meshFilter.sharedMesh.GetOctree();
    }

    private void OnDrawGizmos()
    {
        if (meshFilter == null) return;
        if (meshFilter.sharedMesh == null) return;
        if (meshOctree == null)
        {
            meshOctree = meshFilter.sharedMesh.GetOctree();
        }

        Gizmos.color = Color.yellow;
        Gizmos.DrawSphere(transform.position, 0.1f);

        Vector3 origin = transform.position;
        Vector3 dir = transform.forward;

        // Convert to local space
        var matrix = meshFilter.transform.worldToLocalMatrix;

        origin = matrix * origin;
        dir = matrix * new Vector4(dir.x, dir.y, dir.z, 0.0f);

        Triangle triangle = null;
        float t = float.MaxValue;
        if (meshOctree.Raycast(origin, dir, maxDistance, ref triangle, ref t))
        {
            Gizmos.color = Color.red;
            Gizmos.DrawLine(transform.position, transform.position + transform.forward * t);

            matrix = meshFilter.transform.localToWorldMatrix;

            Gizmos.color = Color.magenta;
            Gizmos.DrawLine(matrix * triangle.GetVertex(0), matrix * triangle.GetVertex(1));
            Gizmos.DrawLine(matrix * triangle.GetVertex(1), matrix * triangle.GetVertex(2));
            Gizmos.DrawLine(matrix * triangle.GetVertex(2), matrix * triangle.GetVertex(0));
        }
        else
        {
            Gizmos.color = Color.green;
            Gizmos.DrawLine(transform.position, transform.position + transform.forward * maxDistance);

            if (meshFilter.sharedMesh.Raycast(origin, dir, maxDistance, out submeshId, out triId))
            {
                triangle = meshFilter.sharedMesh.GetTriangle(submeshId, triId);

                Gizmos.color = Color.cyan;
                Gizmos.DrawLine(matrix * triangle.GetVertex(0), matrix * triangle.GetVertex(1));
                Gizmos.DrawLine(matrix * triangle.GetVertex(1), matrix * triangle.GetVertex(2));
                Gizmos.DrawLine(matrix * triangle.GetVertex(2), matrix * triangle.GetVertex(0));

            }
        }

        Gizmos.color = new Color(1.0f, 1.0f, 0.0f, 0.5f);
        meshOctree.DrawGizmos(50);
    }
}
