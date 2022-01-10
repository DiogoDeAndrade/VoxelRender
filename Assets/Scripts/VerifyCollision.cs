using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class VerifyCollision : MonoBehaviour
{
    public Vector3      size = Vector3.one;
    public MeshFilter   targetMesh;
    public VoxelObject  targetVO;

    public Vector3Int   voxelCoordinates;
    public int          voxelValue;
    public int          collisionSubmesh = -1;
    public int          collisionTriangle = -1;
    public Vector3      voxelAABBMin;
    public Vector3      voxelAABBMax;

    // Start is called before the first frame update
    void Start()
    {
        
    }

    // Update is called once per frame
    void Update()
    {
        
    }

    private void OnDrawGizmosSelected()
    {
        Vector3 aabb_min = Vector3.zero;
        Vector3 aabb_max = Vector3.zero;

        if (targetVO)
        {
            voxelCoordinates = new Vector3Int(-1, -1, -1);

            Vector3 center = transform.position - targetMesh.transform.position - targetVO.offset;

            voxelCoordinates.x = Mathf.FloorToInt(center.x / targetVO.voxelSize.x);
            voxelCoordinates.y = Mathf.FloorToInt(center.y / targetVO.voxelSize.y);
            voxelCoordinates.z = Mathf.FloorToInt(center.z / targetVO.voxelSize.z);

            (aabb_min, aabb_max) = targetVO.GetAABB(voxelCoordinates);

            Gizmos.color = new Color(1.0f, 1.0f, 0.0f, 0.5f);
            Gizmos.DrawCube((aabb_min + aabb_max) * 0.5f, (aabb_max - aabb_min));

            if (targetVO.data.Length == 0)
            {
                Debug.LogWarning("No data on voxel object!");
            }
            else if ((voxelCoordinates.x >= 0) && (voxelCoordinates.x < targetVO.gridSize.x) &&
                     (voxelCoordinates.y >= 0) && (voxelCoordinates.y < targetVO.gridSize.y) &&
                     (voxelCoordinates.z >= 0) && (voxelCoordinates.z < targetVO.gridSize.z))
            {
                voxelValue = targetVO.Get(voxelCoordinates.x, voxelCoordinates.y, voxelCoordinates.z);
            }
            else
            {
                voxelValue = -1;
            }
        }

        if (targetMesh == null)
        {
            Gizmos.color = Color.magenta;
        }
        else
        {
            Gizmos.color = Color.green;
            collisionTriangle = -1;
            collisionSubmesh = -1;

            Mesh mesh = targetMesh.sharedMesh;
            if (mesh)
            {
                aabb_min = transform.position - targetMesh.transform.position - size * 0.5f;
                aabb_max = transform.position - targetMesh.transform.position + size * 0.5f;
                voxelAABBMin = aabb_min;
                voxelAABBMax = aabb_max;
                var vertex = mesh.vertices;
                for (int submesh = 0; submesh < mesh.subMeshCount; submesh++)
                {
                    var indices = mesh.GetTriangles(submesh);

                    for (int i = 0; i < indices.Length; i+=3)
                    {
                        var v1 = vertex[indices[i]];
                        var v2 = vertex[indices[i+1]];
                        var v3 = vertex[indices[i+2]];

                        if (AABB.Intersects(aabb_min, aabb_max, v1, v2, v3))
                        {
                            Gizmos.color = Color.red;
                            collisionTriangle = i;
                            collisionSubmesh = submesh;
                            break;
                        }
                    }

                    if (collisionTriangle != -1) break;
                }
            }
        }
        Gizmos.DrawWireCube(transform.position, size);
    }
}
