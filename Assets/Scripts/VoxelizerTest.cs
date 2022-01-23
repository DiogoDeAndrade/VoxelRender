using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using NaughtyAttributes;

public class VoxelizerTest : MonoBehaviour
{
    public MeshFilter       sourceMesh;
    public VoxelObject      voxelObject;
    public MeshFilter       navMeshDisplay;
    public BoundaryDisplay  boundaryDisplay;
    public BoundaryDisplay  simplifiedBoundaryDisplay;
    public bool             displayGrid;
    [Range(0.1f, 40.0f)]
    public float        density = 1.0f;
    public float        triangleScale = 1.0f;
    public float        gridScale = 1.0f;
    public bool         markWalkable = false;
    [ShowIf("markWalkable")]
    public float        agentHeight = 2.0f;
    [ShowIf("markWalkable")]
    public float        agentStep = 0.2f;
    [ShowIf("markWalkable")]
    public float        minRegionAreaPercentage = 0.1f;
    [ShowIf("markWalkable")]
    public bool         buildNavMesh = true;
    [ShowIf(EConditionOperator.And, "markWalkable", "buildNavMesh")]
    public bool         simplifyNavMesh = true;
    [ShowIf(EConditionOperator.And, "markWalkable", "buildNavMesh")]
    public bool         simplifyBoundary = true;
    [ShowIf(EConditionOperator.And, "markWalkable", "simplifyBoundary")]
    public float        boundarySimplificationMaxDistance = 0.0f;

    public bool         displayVertex = false;
    [ShowIf("displayVertex")]
    public int          vertexId = 0;

    [Button("Voxelize")]
    void Voxelize()
    {
        System.Diagnostics.Stopwatch stopwatch = System.Diagnostics.Stopwatch.StartNew();

        Mesh mesh = sourceMesh.sharedMesh;

        var voxelData = VoxelTools.Voxelize(mesh, density, triangleScale, gridScale);

        var t0 = stopwatch.ElapsedMilliseconds;
        Debug.Log("Voxelization time = " + t0);

        boundaryDisplay?.Clear();

        if (markWalkable)
        {
            int voxelStep = Mathf.FloorToInt(agentStep / voxelData.voxelSize.y);
            int voxelMaxHeight = Mathf.CeilToInt(agentHeight / voxelData.voxelSize.y);

            Debug.Log("Voxel Height = " + voxelMaxHeight);
            Debug.Log("Voxel Step = " + voxelStep);

            t0 = stopwatch.ElapsedMilliseconds;
            VoxelTools.MarkMinHeight(voxelData, voxelMaxHeight, 2);
            Debug.Log("Mark min height = " + (stopwatch.ElapsedMilliseconds - t0));

            int countWalkable = VoxelTools.CountVoxel(voxelData, 2);
            Debug.Log("Walkable area = " + countWalkable + " voxels");

            t0 = stopwatch.ElapsedMilliseconds;
            var regions = VoxelTools.MarkRegions(voxelData, voxelStep, 2, 32, 36);
            Debug.Log("Mark regions = " + (stopwatch.ElapsedMilliseconds - t0));

            if (regions != null)
            {
                List<int> validVoxels = new List<int>();

                t0 = stopwatch.ElapsedMilliseconds;
                for (int i = 0; i < regions.Count; i++)
                {
                    var r = regions[i];
                    float p = r.size / (float)countWalkable;
                    
                    if (p < minRegionAreaPercentage)
                    {
                        Debug.Log("(X) Region " + i + " : Size = " + r.size + " = " + (p * 100.0f) + "%");

                        VoxelTools.FloodFillWithStep(voxelData, voxelStep, r.startX, r.startY, r.startZ, r.voxelId, 1);
                    }
                    else
                    {
                        Debug.Log("( ) Region " + i + " : Size = " + regions[i].size + " = " + (p * 100.0f) + "%");

                        validVoxels.Add(r.voxelId);
                    }
                }
                Debug.Log("Remove small regions = " + (stopwatch.ElapsedMilliseconds - t0));

                if (buildNavMesh)
                {
                    t0 = stopwatch.ElapsedMilliseconds;

                    VoxelNavMesh vnm = new VoxelNavMesh();

                    vnm.Build(voxelData, validVoxels,
                              voxelStep * voxelData.voxelSize.y,
                              simplifyNavMesh, (simplifyNavMesh)?(true):(false));

                    if (navMeshDisplay)
                    {
                        navMeshDisplay.sharedMesh = vnm.GetMesh();
                    }

                    Debug.Log("Build navmesh = " + (stopwatch.ElapsedMilliseconds - t0));

                    if (boundaryDisplay)
                    {
                        var unsimplifiedMesh = (simplifyNavMesh)?(vnm.GetUnsimplifiedMesh()):(vnm.GetMesh());
                        if (unsimplifiedMesh)
                        {
                            Topology topology = new Topology(unsimplifiedMesh);

                            boundaryDisplay.boundary = topology.GetBoundary();

                            if (simplifyBoundary)
                            {
                                if (simplifiedBoundaryDisplay)
                                {
                                    simplifiedBoundaryDisplay.boundary = topology.GetBoundary();
                                    simplifiedBoundaryDisplay.boundary.Simplify(boundarySimplificationMaxDistance);
                                }
                                else
                                {
                                    boundaryDisplay.boundary.Simplify(boundarySimplificationMaxDistance);
                                }
                            }
                        }
                    }
                }
            }
        }

        t0 = stopwatch.ElapsedMilliseconds;

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
        Debug.Log("Mesh generation time = " + (stopwatch.ElapsedMilliseconds - t0));
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

        if (displayVertex)
        {
            var prevMatrix = Gizmos.matrix;
            Gizmos.matrix = navMeshDisplay.transform.localToWorldMatrix;

            Mesh mesh = navMeshDisplay.sharedMesh;
            var  vertices = mesh.vertices;

            Gizmos.color = Color.yellow;
            Gizmos.DrawSphere(vertices[vertexId], 0.02f);

            Gizmos.matrix = prevMatrix;
        }
    }
}
