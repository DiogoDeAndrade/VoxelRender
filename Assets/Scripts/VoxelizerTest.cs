using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using NaughtyAttributes;

public class VoxelizerTest : MonoBehaviour
{
    public MeshFilter   sourceMesh;
    public VoxelObject  voxelObject;
    public bool         displayGrid;
    [Range(0.1f, 40.0f)]
    public float        density = 1.0f;
    public float        triangleScale = 1.0f;
    public float        gridScale = 1.0f;

    [Button("Voxelize")]
    void Voxelize()
    {
        System.Diagnostics.Stopwatch stopwatch = System.Diagnostics.Stopwatch.StartNew();

        Mesh mesh = sourceMesh.sharedMesh;

        var voxelData = MeshTools.Voxelize(mesh, density, triangleScale, gridScale);

        voxelObject.gridSize = voxelData.gridSize;
        voxelObject.voxelSize = voxelData.voxelSize;
        voxelObject.offset = voxelData.offset;
        voxelObject.data = voxelData.data;

        MeshFilter targetMeshFilter = voxelObject.GetComponent<MeshFilter>();
        if (targetMeshFilter)
        {
            targetMeshFilter.mesh = voxelObject.GetMesh();
        }

        stopwatch.Stop();
        Debug.Log("Voxelization time = " + stopwatch.ElapsedMilliseconds);
    }

    private void OnDrawGizmosSelected()
    {
        if (!sourceMesh) return;

        if (displayGrid)
        {
            var prevMatrix = Gizmos.matrix;
            Gizmos.matrix = sourceMesh.transform.localToWorldMatrix;

            Mesh   mesh = sourceMesh.sharedMesh;
            Bounds bounds = mesh.bounds;
            bounds.Expand(bounds.size * gridScale);

            Vector3Int gridSize = new Vector3Int(Mathf.CeilToInt(bounds.size.x * density),
                                                 Mathf.CeilToInt(bounds.size.y * density),
                                                 Mathf.CeilToInt(bounds.size.z * density));
            float voxelSize = 1.0f / density;

            for (int z = 0; z < gridSize.z; z++)
            {
                for (int y = 0; y < gridSize.y; y++)
                {
                    for (int x = 0; x < gridSize.x; x++)
                    {
                        var center = new Vector3(x * voxelSize + bounds.min.x + voxelSize * 0.5f, y * voxelSize + bounds.min.y + voxelSize * 0.5f, z * voxelSize + bounds.min.z + voxelSize * 0.5f);
                        Gizmos.DrawWireCube(center, new Vector3(voxelSize, voxelSize, voxelSize));
                    }
                }
            }

            Gizmos.matrix = prevMatrix;
        }
    }
}
