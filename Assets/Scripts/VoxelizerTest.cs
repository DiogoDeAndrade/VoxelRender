using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using NaughtyAttributes;
using static RcdtcsUnityUtils;

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
    public bool         useRecast = false;
    private bool        useCustom => !useRecast;
    [ShowIf(EConditionOperator.And, "markWalkable", "buildNavMesh", "useCustom")]
    public bool         simplifyNavMesh = true;
    [ShowIf(EConditionOperator.And, "markWalkable", "buildNavMesh", "useCustom")]
    public VoxelNavMesh.TriangulationMesh   simplificationMethod;
    [ShowIf(EConditionOperator.And, "markWalkable", "buildNavMesh", "useCustom")]
    public bool         simplifyBoundary = true;
    [ShowIf(EConditionOperator.And, "markWalkable", "buildNavMesh", "simplifyBoundary", "useCustom")]
    public float        boundarySimplificationMaxDistance = 0.0f;
    [ShowIf(EConditionOperator.And, "markWalkable", "buildNavMesh", "useRecast")]
    public float        minArea = 1.0f;
    [ShowIf(EConditionOperator.And, "markWalkable", "buildNavMesh", "useRecast")]
    public float        agentRadius = 0.1f;
    [ShowIf(EConditionOperator.And, "markWalkable", "buildNavMesh", "useRecast")]
    public bool         detectInOut = false;
    [ShowIf(EConditionOperator.And, "markWalkable", "buildNavMesh", "useRecast", "detectInOut")]
    public float        inOutRange = 0.1f;
    [ShowIf(EConditionOperator.And, "markWalkable", "buildNavMesh", "useRecast", "detectInOut")]
    public Vector3      inOutOffsetScale = new Vector3(0.0f, 0.5f, 0.0f);
    [ShowIf(EConditionOperator.And, "markWalkable", "buildNavMesh", "useRecast", "detectInOut")]
    public int          inOutProbeCount = 16;
    [ShowIf(EConditionOperator.And, "markWalkable", "buildNavMesh", "useRecast", "detectInOut")]
    public bool         inOutProbeInsideCheck = true;
    [ShowIf(EConditionOperator.And, "markWalkable", "buildNavMesh", "useRecast", "detectInOut")]
    public bool         detectInOutProbeDebugView = false;
    [ShowIf(EConditionOperator.And, "markWalkable", "buildNavMesh", "useRecast", "detectInOut", "detectInOutProbeDebugView")]
    public int          detectInOutProbeDebugViewEdgeIndex = -1;
    [ShowIf(EConditionOperator.And, "markWalkable", "buildNavMesh", "useRecast", "detectInOut")]
    public bool         displayInOut = true;

    public bool         displayVertex = false;
    [ShowIf("displayVertex")]
    public int          vertexId = 0;  

    [System.Serializable]
    struct Edge
    {
        public Vector3 p0;
        public Vector3 p1;
    }

    [SerializeField, HideInInspector] private List<Edge> inOutEdge;
     
    [Button("Build")]
    void Voxelize()
    {
        System.Diagnostics.Stopwatch stopwatch = System.Diagnostics.Stopwatch.StartNew();

        long t0 = 0;

        if (useRecast)
        {
            voxelObject.data = null;

            t0 = stopwatch.ElapsedMilliseconds;

            if (markWalkable && buildNavMesh)
            {
                RecastMeshParams navMeshParams = new RecastMeshParams();
                navMeshParams.m_agentHeight = agentHeight;
                navMeshParams.m_agentRadius = agentRadius;
                navMeshParams.m_agentMaxClimb = agentStep;
                navMeshParams.m_agentMaxSlope = 45;

                navMeshParams.m_cellSize = 1.0f / density;
                navMeshParams.m_cellHeight = 1.0f / density;

                navMeshParams.m_regionMinSize = minArea;
                navMeshParams.m_regionMergeSize = 0;
                navMeshParams.m_monotonePartitioning = false;

                navMeshParams.m_edgeMaxLen = 1;
                navMeshParams.m_edgeMaxError = 1;
                navMeshParams.m_vertsPerPoly = 6;
                navMeshParams.m_detailSampleDist = 1;
                navMeshParams.m_detailSampleMaxError = 1;

                SystemHelper recast = new SystemHelper();

                recast.SetNavMeshParams(navMeshParams);
                recast.ClearComputedData();
                recast.ClearMesh();
                recast.AddMesh(sourceMesh.sharedMesh, sourceMesh.gameObject);
                recast.ComputeSystem();

                var t1 = stopwatch.ElapsedMilliseconds;
                Debug.Log("Navmesh generation = " + t1);

                Mesh navMesh = recast.GetNavMesh(transform.worldToLocalMatrix);
                if (navMeshDisplay)
                {
                    navMeshDisplay.sharedMesh = navMesh;
                }

                if (simplifiedBoundaryDisplay)
                {
                    simplifiedBoundaryDisplay.boundary = null;
                }
                if (boundaryDisplay)
                {
                    boundaryDisplay.boundary = null;

                    Topology topology = new Topology(navMesh);

                    boundaryDisplay.boundary = topology.GetBoundary();
                }

                if (detectInOut)
                {
                    DetectInOut();
                }
            }
        }
        else
        {
            Mesh mesh = sourceMesh.sharedMesh;

            var voxelData = VoxelTools.Voxelize(mesh, density, triangleScale, gridScale);

            t0 = stopwatch.ElapsedMilliseconds;
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
                                  simplificationMethod, simplifyBoundary,
                                  voxelStep * voxelData.voxelSize.y,
                                  simplifyNavMesh, (simplifyNavMesh) ? (true) : (false));

                        if (navMeshDisplay)
                        {
                            navMeshDisplay.sharedMesh = vnm.GetMesh();
                        }

                        Debug.Log("Build navmesh = " + (stopwatch.ElapsedMilliseconds - t0));

                        if (boundaryDisplay)
                        {
                            var unsimplifiedMesh = (simplifyNavMesh) ? (vnm.GetUnsimplifiedMesh()) : (vnm.GetMesh());
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
        }

        MeshFilter targetMeshFilter = voxelObject.GetComponent<MeshFilter>();
        if (targetMeshFilter)
        {
            targetMeshFilter.mesh = voxelObject.GetMesh();
        }

        stopwatch.Stop();
        Debug.Log("Mesh generation time = " + (stopwatch.ElapsedMilliseconds - t0));
    }

    void DetectInOut()
    {
        inOutEdge = null;

        Boundary boundary;

        if (inOutRange <= 0) return;
        if ((boundaryDisplay == null) || (boundaryDisplay.boundary == null))
        {
            boundary = new Topology(navMeshDisplay.sharedMesh).GetBoundary();
        }
        else
        {
            boundary = boundaryDisplay.boundary;
        }

        if (boundary == null) return;

        inOutEdge = new List<Edge>();

        Mesh mesh = (sourceMesh) ? (sourceMesh.sharedMesh) : (null);

        int triIndex, submeshIndex;

        float step = 2.0f * inOutRange;

        for (int boundaryIndex = 0; boundaryIndex < boundary.Count; boundaryIndex++)
        {
            var polyline = boundary.Get(boundaryIndex);

            for (int edgeCount = 0; edgeCount < polyline.Count; edgeCount++)
            {
                // Get edge
                bool validEdge = true;

                var p1 = polyline[edgeCount];
                var p2 = polyline[(edgeCount + 1) % polyline.Count];

                var edgeDir = p2 - p1;
                var length = edgeDir.magnitude;
                edgeDir = edgeDir / length;

                var count = Mathf.CeilToInt(length / step);
                var inc = length / count;

                var p = p1 + inOutOffsetScale * step;
                for (int i = 0; i < (count + 1); i++)
                {
                    float startAngle = 0;
                    float endAngle = 360;
                    float incAngle = (endAngle - startAngle) / inOutProbeCount;
                    float angle = startAngle;
                    for (int j = 0; j < inOutProbeCount; j++)
                    {
                        Vector3 probeEnd = p + new Vector3(step * Mathf.Cos(Mathf.Deg2Rad * angle),
                                                            0.0f,
                                                            step * Mathf.Sin(Mathf.Deg2Rad * angle));

                        angle += incAngle;

                        if (inOutProbeInsideCheck)
                        {
                            float d1 = Vector3.Dot(probeEnd - p1, edgeDir);
                            float d2 = Vector3.Dot(probeEnd - p2, edgeDir);

                            if ((d1 * d2) > 0)
                            {
                                continue;
                            }
                        }

                        if (mesh)
                        {
                            var d = probeEnd - p;
                            var l = d.magnitude;
                            d /= l;
                            if (mesh.Raycast(p, d, l, out submeshIndex, out triIndex))
                            {
                                validEdge = false;
                                break;
                            }
                        }
                    }

                    if (!validEdge) break;

                    p += edgeDir * inc;
                }

                if (validEdge)
                {
                    inOutEdge.Add(new Edge() { p0 = p1, p1 = p2 });
                }
            }
        }
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

        if ((detectInOut) && (detectInOutProbeDebugView) && (inOutRange > 0))
        {
            Boundary boundary = null;
            Mesh     mesh = (sourceMesh)?(sourceMesh.sharedMesh):(null);

            if ((boundaryDisplay) && (boundaryDisplay.boundary != null))
                boundary = boundaryDisplay.boundary;
            else if ((simplifiedBoundaryDisplay) && (simplifiedBoundaryDisplay.boundary != null))
                boundary = simplifiedBoundaryDisplay.boundary;
            else
            {
                boundary = new Topology(navMeshDisplay.sharedMesh).GetBoundary();
            }
            
            if (boundary != null)
            {
                int triIndex, submeshIndex;
                int edgeId = 0;

                float step = 2.0f * inOutRange;

                var prevMatrix = Gizmos.matrix;
                Gizmos.matrix = sourceMesh.transform.localToWorldMatrix;

                for (int boundaryIndex = 0; boundaryIndex < boundary.Count; boundaryIndex++)
                {
                    var polyline = boundary.Get(boundaryIndex);

                    for (int edgeCount = 0; edgeCount < polyline.Count; edgeCount++)
                    {
                        edgeId++;

                        if (detectInOutProbeDebugViewEdgeIndex != -1)
                        {
                            if (edgeId != detectInOutProbeDebugViewEdgeIndex) continue;
                        }

                        // Get edge
                        var p1 = polyline[edgeCount];
                        var p2 = polyline[(edgeCount + 1) % polyline.Count];

                        var edgeDir = p2 - p1;
                        var length = edgeDir.magnitude;
                        edgeDir = edgeDir / length;

                        var count = Mathf.CeilToInt(length / step);
                        var inc = length / count;

                        var   p = p1 + inOutOffsetScale * step;
                        for (int i = 0; i < (count + 1); i++)
                        {
                            float startAngle = 0;
                            float endAngle = 360;
                            float incAngle = (endAngle - startAngle) / inOutProbeCount;
                            float angle = startAngle;
                            for (int j = 0; j < inOutProbeCount; j++)
                            {
                                Vector3 probeEnd = p + new Vector3(step * Mathf.Cos(Mathf.Deg2Rad * angle),
                                                                   0.0f,
                                                                   step * Mathf.Sin(Mathf.Deg2Rad * angle));

                                angle += incAngle;

                                if (inOutProbeInsideCheck)
                                {
                                    float d1 = Vector3.Dot(probeEnd - p1, edgeDir);
                                    float d2 = Vector3.Dot(probeEnd - p2, edgeDir);

                                    if ((d1 * d2) > 0)
                                    {
                                        continue;
                                    }
                                }

                                if (mesh)
                                {
                                    var d = probeEnd - p;
                                    var l = d.magnitude;
                                    d /= l;
                                    if (mesh.Raycast(p, d, l, out submeshIndex, out triIndex))
                                    {
                                        var triangle = mesh.GetTriangle(submeshIndex, triIndex);

                                        Gizmos.color = Color.yellow;
                                        Gizmos.DrawLine(triangle.GetVertex(0), triangle.GetVertex(1));
                                        Gizmos.DrawLine(triangle.GetVertex(1), triangle.GetVertex(2));
                                        Gizmos.DrawLine(triangle.GetVertex(2), triangle.GetVertex(0));

                                        Gizmos.color = Color.red;
                                    }
                                    else Gizmos.color = Color.green;
                                }
                                else Gizmos.color = Color.yellow;

                                Gizmos.DrawLine(p, probeEnd);
                            }

                            p += edgeDir * inc;
                        }

                        edgeId++;
                    }
                }

                Gizmos.matrix = prevMatrix;
            }
        }
        
        if ((displayInOut) && (inOutEdge != null))
        {
            var prevMatrix = UnityEditor.Handles.matrix;
            UnityEditor.Handles.matrix = navMeshDisplay.transform.localToWorldMatrix;

            foreach (var edge in inOutEdge)
            {
                UnityEditor.Handles.DrawBezier(edge.p0, edge.p1, edge.p0, edge.p1, Color.red, null, 15.0f);
            }

            UnityEditor.Handles.matrix = prevMatrix;
        }
    }
}
