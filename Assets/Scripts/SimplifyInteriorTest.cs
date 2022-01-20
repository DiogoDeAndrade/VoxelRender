using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using NaughtyAttributes;

[RequireComponent(typeof(MeshFilter))]
[RequireComponent(typeof(MeshRenderer))]
public class SimplifyInteriorTest : MonoBehaviour
{
    public enum MeshType { Test1, Test2, QuadGrid, QuadGridX2, QuadTriMinimal };

    public MeshType type;

    [Button("BuildMesh")]
    void BuildMesh()
    {
        Mesh mesh = null;
        if (type == MeshType.Test1)
        {
            mesh = new Mesh();
            mesh.name = "Test1";

            // Test 1 - Triangle with extra vertice inside (3 tris -> 1 tri)
            mesh.vertices = new Vector3[] {
                new Vector3(0.0f, 0.0f, 0.0f),
                new Vector3(1.0f, 0.0f, 1.0f),
                new Vector3(2.0f, 0.0f, 0.0f),
                new Vector3(1.0f, 0.0f, 2.0f),
            };
            mesh.triangles = new int[]
            {
                0, 1, 2,
                1, 3, 2,
                0, 3, 1
            };
        }
        else if (type == MeshType.Test2)
        {
            mesh = new Mesh();
            mesh.name = "Test2";

            // Test 2 - 4 quads in a grid
            mesh.vertices = new Vector3[] {
                new Vector3(0.0f, 0.0f, 0.0f),
                new Vector3(1.0f, 0.0f, 0.0f),
                new Vector3(2.0f, 0.0f, 0.0f),
                new Vector3(0.0f, 0.0f, 1.0f),
                new Vector3(1.0f, 0.0f, 1.0f),
                new Vector3(2.0f, 0.0f, 1.0f),
                new Vector3(0.0f, 0.0f, 2.0f),
                new Vector3(1.0f, 0.0f, 2.0f),
                new Vector3(2.0f, 0.0f, 2.0f),
                new Vector3(-1.0f, 0.0f, 1.0f),
                new Vector3(-1.0f, 0.0f, 2.0f),
            };
            mesh.triangles = new int[]
            {
                0, 3, 4, 0, 4, 1,
                1, 4, 5, 1, 5, 2,
                3, 6, 7, 3, 7, 4,
                4, 7, 8, 4, 8, 5,
                9, 10, 6, 9, 6, 3
            };
        }
        else if (type == MeshType.QuadGrid)
        {
            mesh = GeometricFactory.BuildPlane("QuadGrid", 2.0f, 20, Vector2.one, false);
        }
        else if (type == MeshType.QuadGridX2)
        {
            mesh = GeometricFactory.BuildPlane("QuadGrid", 2.0f, 40, Vector2.one, false);
        }
        else if (type == MeshType.QuadTriMinimal)
        {
            mesh = GeometricFactory.BuildPlane("QuadGrid", 2.0f, 2, Vector2.one, false);
        }

        GetComponent<MeshFilter>().sharedMesh = mesh;
    }

    [Button("Simplify")]
    void SimplifyMesh()
    {
        System.Diagnostics.Stopwatch stopwatch = System.Diagnostics.Stopwatch.StartNew();

        var t0 = stopwatch.ElapsedMilliseconds;
        var meshFilter = GetComponent<MeshFilter>();
        meshFilter.sharedMesh = MeshTools.SimplifyMeshInterior(meshFilter.sharedMesh, 0.001f);
        Debug.Log("Simplify time = " + (stopwatch.ElapsedMilliseconds - t0));
    }


    [Button("Manual Simplify")]
    void ManualSimplifyMesh()
    {
        var meshFilter = GetComponent<MeshFilter>();
        var mesh = meshFilter.sharedMesh;

        var topology = new Topology(mesh);

        topology.CollapseEdge(4, 3);

        meshFilter.sharedMesh = MeshTools.FromTopology(topology);
    }
}
