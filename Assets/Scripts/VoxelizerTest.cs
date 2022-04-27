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
    [Range(1.0f, 5.0f)]
    public float        verticalDensityScale = 1.0f;
    public bool         markWalkable = false;
    [ShowIf("markWalkable")]
    public float        agentRadius = 0.1f;
    [ShowIf("markWalkable")]
    public float        agentHeight = 2.0f;
    [ShowIf("markWalkable")]
    public float        agentStep = 0.2f;
    [ShowIf("markWalkable")]
    public float        agentMaxSlope = 45.0f;
    [ShowIf("markWalkable")]
    public bool         buildNavMesh = true;
    [ShowIf(EConditionOperator.And, "markWalkable", "buildNavMesh")]
    public bool         useRecast = false;
    private bool        useCustom => !useRecast && markWalkable && buildNavMesh;
    [ShowIf("useCustom")]
    public float        triangleScale = 1.0f;
    [ShowIf("useCustom")]
    public float        gridScale = 1.0f;
    [ShowIf("useCustom")]
    public float        minRegionAreaPercentage = 0.1f;
    [ShowIf("useCustom")]
    public bool         simplifyNavMesh = true;
    [ShowIf("useCustom")]
    public VoxelNavMesh.TriangulationMesh   simplificationMethod;
    [ShowIf("useCustom")]
    public bool         simplifyBoundary = true;
    [ShowIf(EConditionOperator.And, "simplifyBoundary", "useCustom")]
    public float        boundarySimplificationMaxDistance = 0.0f;
    private bool        shouldUseRecast => markWalkable && buildNavMesh && useRecast;
    [ShowIf("shouldUseRecast")]
    public float        minArea = 1.0f;
    [ShowIf("shouldUseRecast")]
    public bool         detectInOut = false;
    private bool        shouldDetectInOut => markWalkable && buildNavMesh && useRecast && detectInOut;
    [ShowIf("shouldDetectInOut")]
    public bool         inOutPerEdge = true;
    private bool        inOutContinuous => !inOutPerEdge && detectInOut && markWalkable && buildNavMesh && useRecast;
    [ShowIf("inOutContinuous")]
    public float        inOutSampleLength = 0.01f;
    [ShowIf("inOutContinuous")]
    public float        inOutSampleDirTolerance = 5;
    [ShowIf("inOutContinuous")]
    public bool         raycastParallelEdgePlane = true;
    [ShowIf("inOutContinuous")]
    public float        edgeMaxAngleWithXZ = 180;
    [ShowIf("shouldDetectInOut")]
    public float        inOutRange = 0.1f;
    [ShowIf("shouldDetectInOut")]
    public int          inOutProbeCount = 16;
    [ShowIf("shouldDetectInOut"), Range(0.0f, 1.0f)]
    public float        probeOffsetScale = 0.1f;
    [ShowIf("shouldDetectInOut")]
    public float        probeRadiusScale = 1.5f;
    [ShowIf("shouldDetectInOut")]
    public bool         inOutProbeInsideCheck = true;
    [ShowIf("shouldDetectInOut")]
    public bool         removeInteriorEdges = true;
    [ShowIf(EConditionOperator.And, "shouldDetectInOut", "removeInteriorEdges")]
    public bool         perpendicularTest = true;
    private bool        circularTest => !perpendicularTest && removeInteriorEdges && shouldDetectInOut;
    [ShowIf("circularTest"), Range(1.0f, 4.0f)]
    public float        raycastInteriorEdgeExcentricity = 1.0f;
    [ShowIf("shouldDetectInOut")]
    public bool         detectInOutProbeDebugView = false;
    [ShowIf(EConditionOperator.And, "shouldDetectInOut", "detectInOutProbeDebugView")]
    public int          detectInOutProbeDebugViewEdgeIndex = -1;
    [ShowIf("shouldDetectInOut")]
    public bool         displayInOut = true;
    [ShowIf("shouldDetectInOut")]
    public bool         displayInOutNormals = false;
    [ShowIf("shouldDetectInOut")]
    public bool         matchSourceGeometry = true;
    private bool        shouldMatchSourceGeometry => matchSourceGeometry && shouldDetectInOut;
    [ShowIf("shouldMatchSourceGeometry")]
    public float        directionAngleTolerance = 5;
    [ShowIf("shouldMatchSourceGeometry")]
    public float        normalAngleTolerance = 45;

    public bool         displayVertex = false;
    [ShowIf("displayVertex")]
    public int          vertexId = 0;  

    public bool         debugRayEnabled = false;
    [ShowIf(EConditionOperator.And, "debugRayEnabled", "removeInteriorEdges")]
    public bool         debugRayExteriorTest = false;
    
    public bool         debugSamplePoints = false;
    [ShowIf("debugSamplePoints")]
    public int          debugSampleId = -1;

    public bool         debugEdges = false;
    [ShowIf("debugEdges")]
    public int          debugEdgeId = -1;
    [ShowIf("debugEdges")]
    public int          debugSubEdgeId = -1;

    private List<Vector3>   samplePoints;



    [System.Serializable]
    struct Edge
    {
        public Vector3 p0;
        public Vector3 p1;
        public Vector3 normal;
        public Vector3 perp;
    }

    [SerializeField, HideInInspector] private List<Edge> inOutEdge;

    struct DebugRay
    {
        public  Color   color;
        public  Vector3 start;
        public  Vector3 dir;
        public  float   dist;
        public  bool    hit;
        public  float   t;
        public  Mesh    mesh;
        public  int     submeshId;
        public  int     triId;
    }
    private List<DebugRay> debugRays = new List<DebugRay>();

    struct DebugEdge
    {
        public Color    color;
        public Vector3  start;
        public Vector3  end;
        public int      edgeId;
        public int      subEdgeId;
    }
    private List<DebugEdge> debugEdgeList = new List<DebugEdge>();

    [Button("Build")]
    void Voxelize()
    {
        debugRays.Clear();

        System.Diagnostics.Stopwatch stopwatch = System.Diagnostics.Stopwatch.StartNew();

        long t0 = 0;

        if (useRecast)
        {
            voxelObject.data = null;

            t0 = stopwatch.ElapsedMilliseconds;

            if (markWalkable && buildNavMesh && (sourceMesh != null))
            {
                UnityEditor.EditorUtility.DisplayProgressBar("Building...", "Creating nav mesh...", 0.0f);

                RecastMeshParams navMeshParams = new RecastMeshParams();
                navMeshParams.m_cellSize = 1.0f / density;
                navMeshParams.m_cellHeight = 1.0f / (density * verticalDensityScale);

                navMeshParams.m_agentHeight = agentHeight;
                navMeshParams.m_agentRadius = agentRadius;
                navMeshParams.m_agentMaxClimb = agentStep;
                navMeshParams.m_agentMaxSlope = agentMaxSlope;

                Debug.Log($"Voxel size = {navMeshParams.m_cellSize}, {navMeshParams.m_cellHeight}, {navMeshParams.m_cellSize}");

                navMeshParams.m_regionMinSize = minArea;
                navMeshParams.m_regionMergeSize = minArea * 10.0f;
                navMeshParams.m_monotonePartitioning = false;

                navMeshParams.m_edgeMaxLen = 0.5f;
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

                Mesh navMesh = recast.GetPolyMesh(transform.worldToLocalMatrix);
                if (navMeshDisplay)
                {
                    navMeshDisplay.sharedMesh = navMesh;
                }

                UnityEditor.EditorUtility.DisplayProgressBar("Building...", "Bulding boundaries...", 0.2f);

                if (simplifiedBoundaryDisplay)
                {
                    simplifiedBoundaryDisplay.boundary = null;
                }
                if (boundaryDisplay)
                {
                    boundaryDisplay.boundary = null;

                    Topology topology = new Topology(navMesh);
                    topology.ComputeTriangleNormals();

                    boundaryDisplay.boundary = topology.GetBoundary();                    
                }

                if (detectInOut)
                {
                    UnityEditor.EditorUtility.DisplayProgressBar("Building...", "Detecting in/out...", 0.25f);

                    if (inOutPerEdge)
                    {
                        DetectInOutPerEdge();
                    }
                    else
                    {
                        DetectInOutContinuous();
                    }

                    UnityEditor.EditorUtility.DisplayProgressBar("Building...", "Matching source geometry...", 0.75f);

                    if (matchSourceGeometry)
                    {
                        MatchSourceGeometry();
                    }
                }
                else
                {
                    inOutEdge.Clear();
                }

                UnityEditor.EditorUtility.ClearProgressBar();
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
                                topology.ComputeTriangleNormals();

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

    void DetectInOutPerEdge()
    {
        inOutEdge = null;

        Boundary boundary;

        if (inOutRange <= 0) return;
        if ((boundaryDisplay == null) || (boundaryDisplay.boundary == null))
        {
            var topology = new Topology(navMeshDisplay.sharedMesh);
            topology.ComputeTriangleNormals();
            boundary = topology.GetBoundary();
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

        int totalEdges = 0;
        for (int boundaryIndex = 0; boundaryIndex < boundary.Count; boundaryIndex++)
        {
            var polyline = boundary.Get(boundaryIndex);
            totalEdges += polyline.Count;
        }
    
        int edgeId = 0;

        for (int boundaryIndex = 0; boundaryIndex < boundary.Count; boundaryIndex++)
        {
            var polyline = boundary.Get(boundaryIndex);

            for (int edgeCount = 0; edgeCount < polyline.Count; edgeCount++)
            {
                edgeId++;

                UnityEditor.EditorUtility.DisplayProgressBar("Building...", "Detecting in/out...", 0.25f + 0.5f * (float)edgeId / totalEdges);

                // Get edge
                bool validEdge = true;

                var p1 = polyline[edgeCount];
                var p2 = polyline[(edgeCount + 1) % polyline.Count];

                var edgeDir = p2 - p1;
                var length = edgeDir.magnitude;
                edgeDir = edgeDir / length;

                var count = Mathf.CeilToInt(length / step);
                var inc = length / count;

                var p = p1 + Vector3.up * (agentStep * probeOffsetScale);
                for (int i = 0; i < (count + 1); i++)
                {
                    float startAngle = 0;
                    float endAngle = 360;
                    float incAngle = (endAngle - startAngle) / inOutProbeCount;
                    float angle = startAngle;
                    for (int j = 0; j < inOutProbeCount; j++)
                    {
                        float   l = Mathf.Max(step, agentRadius * probeRadiusScale);
                        Vector3 probeEnd = p + new Vector3(l * Mathf.Cos(Mathf.Deg2Rad * angle),
                                                           0.0f,
                                                           l* Mathf.Sin(Mathf.Deg2Rad * angle));

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
                            d /= l;
                            bool b = mesh.Raycast(p, d, l, out submeshIndex, out triIndex);

                            if ((detectInOutProbeDebugView) &&
                                ((detectInOutProbeDebugViewEdgeIndex <= 0) || (detectInOutProbeDebugViewEdgeIndex == edgeId)))
                            {
                                debugRays.Add(new DebugRay()
                                {
                                    color = Color.cyan,
                                    start = p,
                                    dir = d,
                                    dist = l,
                                    hit = b,
                                    t = float.MaxValue,
                                    mesh = mesh,
                                    submeshId = submeshIndex,
                                    triId = triIndex
                                });
                            }

                            if (b)
                            {
                                validEdge = false;
                                break;
                            }
                        }
                    }

                    if (!validEdge) break;

                    p += edgeDir * inc;
                }

                if ((validEdge) && (removeInteriorEdges) && (mesh))
                {
                    // Raycast around the center to a certain large distance (10 times the agentRadius for now) and see if at least one ray doesn't intersect anything
                    float startAngle = 0;
                    float endAngle = 360;
                    float incAngle = (endAngle - startAngle) / inOutProbeCount;
                    float angle = startAngle;
                    Vector3 edgeCenter = (p1 + p2) * 0.5f + Vector3.up * (agentStep * probeOffsetScale);

                    validEdge = false;

                    p = edgeCenter;

                    for (int j = 0; j < inOutProbeCount; j++)
                    {
                        float l = agentRadius * 10;
                        Vector3 probeEnd = edgeCenter + new Vector3(l * Mathf.Cos(Mathf.Deg2Rad * angle),
                                                                    0.0f,
                                                                    l * Mathf.Sin(Mathf.Deg2Rad * angle));

                        angle += incAngle;

                        var d = probeEnd - p;
                        d /= l;
                        var b = mesh.Raycast(p, d, l, out submeshIndex, out triIndex);

                        if ((detectInOutProbeDebugView) &&
                            ((detectInOutProbeDebugViewEdgeIndex <= 0) || (detectInOutProbeDebugViewEdgeIndex == edgeId)))
                        {
                            debugRays.Add(new DebugRay()
                            {
                                color = new Color(0.0f, 0.5f, 1.0f),
                                start = p,
                                dir = d,
                                dist = l,
                                hit = b,
                                t = float.MaxValue,
                                mesh = mesh,
                                submeshId = submeshIndex,
                                triId = triIndex
                            });
                        }

                        if (!b)
                        {
                            validEdge = true;
                            break;
                        }

                    }
                }
                
                if (validEdge)
                {
                    var n1 = polyline.GetNormal(edgeCount);
                    var n2 = polyline.GetNormal((edgeCount + 1) % polyline.Count);

                    inOutEdge.Add(new Edge() { p0 = p1, p1 = p2, normal = (n1 + n2) * 0.5f });
                }
            }
        }
    }


    bool CanReachInfiniteCircular(Mesh mesh, Vector3 pt, float range, Vector3 mainDir, Vector3 forwardDir)
    {
        if (mesh == null) return true;

        Vector3 cosDir = Vector3.right;
        Vector3 sinDir = Vector3.forward;

        if (raycastParallelEdgePlane)
        {
            cosDir = mainDir;
            sinDir = forwardDir;
        }

        float startAngle = 0;
        float endAngle = 360;
        float incAngle = (endAngle - startAngle) / inOutProbeCount;
        float angle = startAngle;
        for (int j = 0; j < inOutProbeCount; j++)
        {
            float l = Mathf.Max(range, agentRadius * probeRadiusScale);
            Vector3 probeEnd = pt + cosDir *  l * Mathf.Cos(Mathf.Deg2Rad * angle) + sinDir * l * Mathf.Sin(Mathf.Deg2Rad * angle);

            angle += incAngle;

            float mdp = Vector3.Dot((probeEnd - pt).normalized, mainDir);
            float pdp = Vector3.Dot((probeEnd - pt).normalized, forwardDir);
            if (mdp <= 0.001) continue;
            if (pdp <= 0.001) continue;

            int submeshIndex;
            int triIndex;
            var d = probeEnd - pt;
            d /= l;
            bool b = mesh.Raycast(pt, d, l, out submeshIndex, out triIndex);

            if (debugRayExteriorTest)
            {
                debugRays.Add(new DebugRay()
                {
                    color = Color.blue,
                    start = pt,
                    dir = d,
                    dist = l,
                    hit = b,
                    t = float.MaxValue,
                    mesh = mesh,
                    submeshId = submeshIndex,
                    triId = triIndex
                });//*/
            }

            if (b) return false;
        }

        return true;
    }

    bool CanReachInfinitePerpendicular(Mesh mesh, Vector3 pt0, Vector3 pt1, float range, Vector3 forwardDir)
    {
        if (mesh == null) return true;

        Vector3 pt = pt0;
        Vector3 ptInc = (pt1 - pt0) / inOutProbeCount;

        for (int i = 0; i < inOutProbeCount; i++)
        {
            int submeshIndex;
            int triIndex;

            bool b = mesh.Raycast(pt, forwardDir, range, out submeshIndex, out triIndex);

            if (debugRayExteriorTest)
            {
                debugRays.Add(new DebugRay()
                {
                    color = Color.blue,
                    start = pt,
                    dir = forwardDir,
                    dist = range,
                    hit = b,
                    t = float.MaxValue,
                    mesh = mesh,
                    submeshId = submeshIndex,
                    triId = triIndex
                });//*/
            }

            if (b) return false;

            pt += ptInc;
        }

        return true;
    }

    bool CircularRaycast(Mesh mesh, Vector3 pt, float range, Vector3 edgeDir, Vector3 perpDir, bool isFirst, bool isLast, int sampleId)
    {
        Vector3 cosDir = Vector3.right;
        Vector3 sinDir = Vector3.forward;

        if (raycastParallelEdgePlane)
        {
            cosDir = edgeDir;
            sinDir = perpDir;
        }

        float startAngle = 0;
        float endAngle = 360;
        float incAngle = (endAngle - startAngle) / inOutProbeCount;
        float angle = startAngle;
        for (int j = 0; j < inOutProbeCount; j++)
        {
            float l = Mathf.Max(range, agentRadius * probeRadiusScale);
            Vector3 probeEnd = pt + cosDir * l * Mathf.Cos(Mathf.Deg2Rad * angle) + sinDir * l * Mathf.Sin(Mathf.Deg2Rad * angle);

            angle += incAngle;

            if ((inOutProbeInsideCheck) && (isFirst || isLast) && (edgeDir.magnitude > 0))
            {
                if ((isFirst) && (Vector3.Dot(probeEnd - pt, edgeDir) < 0)) continue;
                if ((isLast) && (Vector3.Dot(probeEnd - pt, edgeDir) > 0)) continue;
            }

            if (mesh)
            {
                int submeshIndex;
                int triIndex;
                var d = probeEnd - pt;
                d /= l;
                bool b = mesh.Raycast(pt, d, l, out submeshIndex, out triIndex);

                if ((debugSamplePoints) &&
                    ((debugSampleId < 0) || (debugSampleId == sampleId)))
                {
                    debugRays.Add(new DebugRay()
                    {
                        color = Color.cyan,
                        start = pt,
                        dir = d,
                        dist = l,
                        hit = b,
                        t = float.MaxValue,
                        mesh = mesh,
                        submeshId = submeshIndex,
                        triId = triIndex
                    });
                }

                if (b) return false;
            }
        }

        return true;
    }

    void DetectInOutContinuous()
    {
        inOutEdge = null;

        Boundary boundary;

        if (inOutRange <= 0) return;
        if ((boundaryDisplay == null) || (boundaryDisplay.boundary == null))
        {
            var topology = new Topology(navMeshDisplay.sharedMesh);
            topology.ComputeTriangleNormals();
            boundary = topology.GetBoundary();
        }
        else
        {
            boundary = boundaryDisplay.boundary;
        }

        if (boundary == null) return;

        inOutEdge = new List<Edge>();

        Mesh mesh = (sourceMesh) ? (sourceMesh.sharedMesh) : (null);

        float step = agentRadius;
        float angleTolerance = Mathf.Cos(inOutSampleDirTolerance * Mathf.Deg2Rad);

        int totalEdges = 0;
        for (int boundaryIndex = 0; boundaryIndex < boundary.Count; boundaryIndex++)
        {
            var polyline = boundary.Get(boundaryIndex);
            totalEdges += polyline.Count;
        }

        int edgeId = 0;
        int sampleId = 0;
        //debugCount = 0;

        samplePoints = new List<Vector3>();

        for (int boundaryIndex = 0; boundaryIndex < boundary.Count; boundaryIndex++)
        {
            var polyline = boundary.Get(boundaryIndex);

            bool isCW = polyline.isCW();

            bool    edgeStarted = false;
            Vector3 candidateEdgeStart = Vector3.zero;
            Vector3 candidateEdgeStartNormal = Vector3.zero;
            int     candidateEdgeStartSampleId = -1;
            Vector3 candidateEdgeEnd = Vector3.zero;
            Vector3 candidateEdgeEndNormal = Vector3.zero;
            int     candidateEdgeEndSampleId = -1;
            Vector3 candidateEdgeDir = Vector3.zero;

            for (int edgeCount = 0; edgeCount < polyline.Count; edgeCount++)
            {
                edgeId++;

                UnityEditor.EditorUtility.DisplayProgressBar("Building...", "Detecting in/out...", 0.25f + 0.5f * (float)edgeId / totalEdges);

                var p1 = polyline[edgeCount];
                var p2 = polyline[(edgeCount + 1) % polyline.Count];
                var n1 = polyline.GetNormal(edgeCount);
                var n2 = polyline.GetNormal((edgeCount + 1) % polyline.Count);

                var edgeDir = p2 - p1;
                var edgeLength = edgeDir.magnitude;
                if (edgeLength == 0) continue;
                edgeDir = edgeDir / edgeLength;

                var nInc = (n2 - n1) / edgeLength;

                var pt = p1;
                var n = n1;
                var nNorm = n.normalized;

                while (Vector3.Dot((pt - p2).normalized, edgeDir) < 0)
                {
                    if ((debugSamplePoints) && ((sampleId == debugSampleId) || (debugSampleId < 0)))
                    {
                        samplePoints.Add(pt);
                    }

                    // Do the raycast for this vertex
                    bool noHit = CircularRaycast(mesh, pt, step, edgeDir, Vector3.Cross(edgeDir, nNorm), !edgeStarted, false, sampleId);
                    if (noHit)
                    {
                        if (!edgeStarted)
                        {
                            edgeStarted = true;
                            candidateEdgeStart = candidateEdgeEnd = pt;
                            candidateEdgeStartNormal = candidateEdgeEndNormal = nNorm;
                            candidateEdgeStartSampleId = sampleId;
                            candidateEdgeDir = Vector3.zero;
                        }
                        else
                        {
                            if (candidateEdgeDir.magnitude == 0)
                            {
                                candidateEdgeEnd = pt;
                                candidateEdgeEndNormal = nNorm;
                                candidateEdgeEndSampleId = sampleId;
                                candidateEdgeDir = (candidateEdgeEnd - candidateEdgeStart).normalized;
                            }
                            else
                            {
                                Vector3 newDir = (pt - candidateEdgeEnd).normalized;

                                if (Vector3.Dot(newDir, candidateEdgeDir) < angleTolerance)
                                {
                                    // Edge direction changed, add the current edge and start a new one
                                    AddInOutEdge(mesh, candidateEdgeStart, candidateEdgeEnd, candidateEdgeStartNormal, candidateEdgeEndNormal, isCW);

                                    candidateEdgeStart = candidateEdgeEnd;
                                    candidateEdgeStartNormal = candidateEdgeEndNormal;
                                    candidateEdgeStartSampleId = candidateEdgeEndSampleId;
                                    candidateEdgeEnd = pt;
                                    candidateEdgeEndNormal = nNorm;
                                    candidateEdgeEndSampleId = sampleId;
                                    candidateEdgeDir = (candidateEdgeEnd - candidateEdgeStart).normalized;
                                }
                                else
                                {
                                    // Continue edge
                                    candidateEdgeEnd = pt;
                                    candidateEdgeEndNormal = nNorm;
                                    candidateEdgeEndSampleId = sampleId;
                                    candidateEdgeDir = (candidateEdgeEnd - candidateEdgeStart).normalized;
                                }
                            }
                        }
                    }
                    else
                    {
                        if (edgeStarted)
                        {
                            // We got a hit, check if it's a potential end
                            var endDir = (pt - candidateEdgeStart);
                            if (endDir.magnitude != 0) endDir.Normalize();

                            if (CircularRaycast(mesh, pt, step, endDir, Vector3.Cross(endDir, nNorm), false, true, sampleId))
                            {
                                // This one is a valid end, add it and reset state to search for a new start
                                candidateEdgeEnd = pt;
                                candidateEdgeEndNormal = nNorm;
                                candidateEdgeEndSampleId = sampleId;
                                candidateEdgeDir = (candidateEdgeEnd - candidateEdgeStart).normalized;

                                AddInOutEdge(mesh, candidateEdgeStart, candidateEdgeEnd, candidateEdgeStartNormal, candidateEdgeEndNormal, isCW);
                            }
                            else
                            {
                                // This one isn't valid as a terminal, so add the current end as the actual end
                                AddInOutEdge(mesh, candidateEdgeStart, candidateEdgeEnd, candidateEdgeStartNormal, candidateEdgeEndNormal, isCW);
                            }

                            // Check if this one can be the start of a new edge
                            if (CircularRaycast(mesh, pt, step, edgeDir, nNorm, true, false, sampleId))
                            {
                                candidateEdgeStart = candidateEdgeEnd = pt;
                                candidateEdgeStartNormal = candidateEdgeEndNormal = nNorm;
                                candidateEdgeStartSampleId = sampleId;
                                candidateEdgeDir = Vector3.zero;
                            }
                            else
                            {
                                edgeStarted = false;
                            }
                        }
                    }

                    pt = pt + edgeDir * inOutSampleLength;
                    n = n + nInc * inOutSampleLength;
                    nNorm = n.normalized;

                    sampleId++;
                }
            }
        }
    }

/*    int debugCount = 0;

    public int addedEdgeTest = 0;
    [ReadOnly] public string rejectionCause;*/

    bool AddInOutEdge(Mesh mesh, Vector3 p0, Vector3 p1, Vector3 n1, Vector3 n2, bool isCW)
    {
        // Check if the edge is larger than the agent radius (smaller size for a connecting edge)
        if (Vector3.Distance(p0, p1) < agentRadius) return false;

        /*debugCount++;
        if ((debugCount != addedEdgeTest) && (addedEdgeTest > 0))
        {
            return false;
        }*/

        //rejectionCause = "No cause";

        Vector3 perp_dir = Vector3.zero;

        if (removeInteriorEdges)
        {
            Vector3 main_dir;
            Vector3 normal = Vector3.up;

            if (raycastParallelEdgePlane)
            {
                normal = ((n1 + n2) * 0.5f).normalized;
            }

            if (perpendicularTest)
            {
                main_dir = (p1 - p0).normalized;
                perp_dir = Vector3.Cross(main_dir, normal);
                //if (isCW) perp_dir = -perp_dir;
                if (!CanReachInfinitePerpendicular(mesh, p0, p1, 4 * agentRadius, perp_dir))
                {
                    //rejectionCause = "Interior segment";
                    return false;
                }
            }
            else
            {
                main_dir = (p1 - p0).normalized;
                perp_dir = Vector3.Cross(main_dir, normal) * raycastInteriorEdgeExcentricity;
                //if (isCW) perp_dir = -perp_dir;
                if (!CanReachInfiniteCircular(mesh, p0, 4 * agentRadius, main_dir, perp_dir))
                {
                    //rejectionCause = "Interior segment (1)";
                    return false;
                }
                main_dir =-main_dir;
                if (!CanReachInfiniteCircular(mesh, p1, 4 * agentRadius, main_dir, perp_dir))
                {
                    //rejectionCause = "Interior segment (2)";
                    return false;
                }
            }
        }

        if (edgeMaxAngleWithXZ < 90)
        {
            Vector3 vec = (p1 - p0).normalized;
            Vector3 projVec = vec.x0z().normalized;

            float dp = Mathf.Abs(Vector3.Dot(projVec, vec));
            if (dp < Mathf.Cos(Mathf.Deg2Rad * edgeMaxAngleWithXZ))
            {
                //rejectionCause = "Angle with XZ plane";
                return false;
            }
        }

        if (perp_dir.magnitude == 0.0f)
        {
            Vector3 normal = Vector3.up;

            if (raycastParallelEdgePlane)
            {
                normal = ((n1 + n2) * 0.5f).normalized;
            }

            var main_dir = (p1 - p0).normalized;
            perp_dir = Vector3.Cross(main_dir, normal);
            if (isCW) perp_dir = -perp_dir;
        }

        inOutEdge.Add(new Edge() { p0 = p0, p1 = p1, normal = (n1 + n2) * 0.5f, perp = perp_dir });

        return true;
    }

    void MatchSourceGeometry()
    {
        Mesh mesh = (sourceMesh) ? (sourceMesh.sharedMesh) : (null);
        if (mesh == null) return;

        var sourceEdges = inOutEdge;
        inOutEdge = new List<Edge>();

        var topology = new Topology(sourceMesh.sharedMesh);
        topology.ComputeTriangleNormals();
        float cosDirectionTolerance = Mathf.Cos(directionAngleTolerance * Mathf.Deg2Rad);
        float cosNormalTolerance = Mathf.Cos(normalAngleTolerance * Mathf.Deg2Rad);
        float toleranceRange = 15.0f * agentRadius;

        int edgeId = 0;

        foreach (var edge in sourceEdges)
        {
            edgeId++;

            if ((debugEdges) && (debugEdgeId != -1) && (debugEdgeId != edgeId)) continue;

            float   score = -float.MaxValue;
            Vector3 candidateStart = Vector3.zero;
            Vector3 candidateEnd = Vector3.zero;
            Vector3 candidateNormal = Vector3.zero;
            Vector3 candidatePerp = Vector3.zero;

            Vector3 edgeStart = edge.p0;
            Vector3 edgeEnd = edge.p1;
            float edgeLength = Vector3.Distance(edgeStart, edgeEnd);
            Vector3 edgeNormal = edge.normal;
            Vector3 edgeCenter = (edgeStart + edgeEnd) * 0.5f;

            debugEdgeList.Add(new DebugEdge() { color = Color.cyan, edgeId = edgeId, subEdgeId = -1, start = edge.p0, end = edge.p1 });

            // Find candidates
            var edgeDir = (edge.p1 - edge.p0).normalized;

            int subEdgeId = 0;

            for (int i = 0; i < topology.nEdges; i++)
            {
                var otherEdge = topology.GetEdge(i);
                var e1 = otherEdge.Item1;
                var e2 = otherEdge.Item2;

                if (Vector3.Distance(e1, edgeStart) > Vector3.Distance(e2, edgeStart))
                {
                    // Switch edges
                    e1 = otherEdge.Item2;
                    e2 = otherEdge.Item1;
                }

                if ((Vector3.Distance(e1, edgeCenter) > toleranceRange) ||
                    (Vector3.Distance(e2, edgeCenter) > toleranceRange)) continue;

                var eC = (e1 + e2) * 0.5f;

                // Check direction tolerance
                float dp = Vector3.Dot((e2 - e1).normalized, edgeDir);
                if ((dp > cosDirectionTolerance) || (-dp > cosDirectionTolerance))
                {
                    // Check normal tolerance
                    dp = Mathf.Abs(Vector3.Dot(topology.GetEdgeNormal(i), edge.normal));

                    if (dp > cosNormalTolerance)
                    {
                        // Check if center of candidate is in front of the current edge
                        if (Vector3.Dot((eC - edgeCenter).normalized, edge.perp) > 0)
                        {
                            subEdgeId++;

                            debugEdgeList.Add(new DebugEdge() { color = Color.green, edgeId = edgeId, subEdgeId = subEdgeId, start = e1, end = e2 });

                            // Project both edges to the XZ plane
                            var edgeStartXZ = edgeStart.x0z();
                            var edgeEndXZ = edgeEnd.x0z();
                            var perp = edge.perp.x0z();
                            var e1XZ = e1.x0z();
                            var e2XZ = e2.x0z();
                            Vector3 intersection;
                            float   t;
                            float   intersectionScore = 0.0f;

                            if (Line.Raycast(edgeStartXZ, perp, toleranceRange, e1XZ, e2XZ, out intersection, out t))
                            {
                                e1 = Line.GetClosestPoint(e1, e2, edgeStart + t * edge.perp);

                                intersectionScore += 4;
                            }
                            if (Line.Raycast(edgeEndXZ, perp, toleranceRange, e1XZ, e2XZ, out intersection, out t))
                            {
                                e2 = Line.GetClosestPoint(e1, e2, edgeEnd + t * edge.perp);

                                intersectionScore += 4;
                            }

                            // The larger the distance, the larger the score
                            float distanceScore = Vector3.Dot(eC.x0z() - edgeStartXZ, edge.perp);
                            // The larger the difference in Y, the lower the score
                            float yScore = -Mathf.Abs(eC.y - edgeCenter.y);
                            float lengthScore = -Mathf.Abs(Vector3.Distance(e1, e2) - edgeLength);

                            float edgeScore = intersectionScore + distanceScore + (yScore * 10) + lengthScore;

                            //Debug.Log($"Score= {edgeScore} / DeltaY = {yScore} / D = {distanceScore} / DiffLength = {lengthScore} / Intersections = {intersectionScore}");

                            if (score < edgeScore)
                            {
                                candidateStart = e1;
                                candidateEnd = e2;
                                candidateNormal = edge.normal;
                                candidatePerp = edge.perp;
                                score = edgeScore;

                                //Debug.Log($"Best score so far: {score}");
                            }

                            debugEdgeList.Add(new DebugEdge() { color = Color.yellow, edgeId = edgeId, subEdgeId = subEdgeId, start = e1, end = e2 });
                        }
                    }
                }
            }

            if (score > -float.MaxValue)
            {
                /*if (raycastEdgeEndpoints)
                {
                    // Do the raycast on both directions along the edge to see if we need to "reduce" the edge to fit
                    Vector3 dir = (candidateEnd - candidateStart).normalized;
                    float candidateLength = (candidateP1 - candidateP0).magnitude;
                    Vector3 candidateCenter = (candidateP1 + candidateP0) * 0.5f;
                    candidateCenter.y += raycastEdgeEndpointsHeight;

                    int submeshId, triHit;
                    float t;
                    bool b1 = mesh.Raycast(candidateCenter, dir, candidateLength * 0.5f, out submeshId, out triHit, out t);

                    if ((detectInOutProbeDebugView) &&
                        ((detectInOutProbeDebugViewEdgeIndex <= 0) || (detectInOutProbeDebugViewEdgeIndex == edgeId)))
                    {
                        debugRays.Add(new DebugRay()
                        {
                            color = new Color(1.0f, 0.5f, 1.0f),
                            start = candidateCenter,
                            dir = dir,
                            dist = candidateLength * 0.5f,
                            hit = b1,
                            t = t,
                            mesh = mesh,
                            submeshId = submeshId,
                            triId = triHit
                        });
                    }

                    bool b2 = mesh.Raycast(candidateCenter, -dir, candidateLength * 0.5f, out submeshId, out triHit, out t);

                    if ((detectInOutProbeDebugView) &&
                        ((detectInOutProbeDebugViewEdgeIndex <= 0) || (detectInOutProbeDebugViewEdgeIndex == edgeId)))
                    {
                        debugRays.Add(new DebugRay()
                        {
                            color = new Color(1.0f, 0.5f, 1.0f),
                            start = candidateCenter,
                            dir = -dir,
                            dist = candidateLength * 0.5f,
                            hit = b1,
                            t = t,
                            mesh = mesh,
                            submeshId = submeshId,
                            triId = triHit
                        });
                    }

                    if (b1)
                    {
                        candidateEnd = candidateCenter + dir * t;
                        candidateEnd.y -= raycastEdgeEndpointsHeight;
                    }

                    if (b2)
                    {
                        candidateStart = candidateCenter - dir * t;
                        candidateStart.y -= raycastEdgeEndpointsHeight;
                    }
                }*/

                inOutEdge.Add(new Edge() { p0 = candidateStart, p1 = candidateEnd, normal = candidateNormal, perp = candidatePerp });
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
                                                 Mathf.CeilToInt(bounds.size.y * density * verticalDensityScale),
                                                 Mathf.CeilToInt(bounds.size.z * density));
            Vector3 voxelSize = new Vector3(1.0f / density, 1.0f / (density * verticalDensityScale), 1.0f / density);

            for (int z = 0; z < gridSize.z; z++)
            {
                for (int y = 0; y < gridSize.y; y++)
                {
                    for (int x = 0; x < gridSize.x; x++)
                    {
                        var center = new Vector3(x * voxelSize.x + bounds.min.x + voxelSize.x * 0.5f, 
                                                 y * voxelSize.y + bounds.min.y + voxelSize.y * 0.5f, 
                                                 z * voxelSize.z + bounds.min.z + voxelSize.z * 0.5f);
                        Gizmos.DrawWireCube(center, new Vector3(voxelSize.x, voxelSize.y, voxelSize.z));
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

        if (debugRayEnabled)
        {
            var prevMatrixHandles = UnityEditor.Handles.matrix;
            var prevMatrixGizmos = Gizmos.matrix;

            UnityEditor.Handles.matrix = navMeshDisplay.transform.localToWorldMatrix;
            Gizmos.matrix = navMeshDisplay.transform.localToWorldMatrix;

            foreach (var ray in debugRays)
            {
                Gizmos.color = ray.color;

                if (ray.hit)
                {
                    Gizmos.color = Color.red;
                }
                Gizmos.DrawLine(ray.start, ray.start + ray.dir * ray.dist);
                if (ray.hit)
                {                    
                    Gizmos.DrawSphere(ray.start + ray.dir * ray.t, 0.05f);

                    if ((ray.mesh) && (ray.submeshId >= 0) && (ray.triId >= 0))
                    {
                        var triangle = ray.mesh.GetTriangle(ray.submeshId, ray.triId);

                        Gizmos.DrawLine(triangle.GetVertex(0), triangle.GetVertex(1));
                        Gizmos.DrawLine(triangle.GetVertex(1), triangle.GetVertex(2));
                        Gizmos.DrawLine(triangle.GetVertex(2), triangle.GetVertex(0));
                    }
                }
            }

            UnityEditor.Handles.matrix = prevMatrixHandles;
            Gizmos.matrix = prevMatrixGizmos;
        }
        
        if ((displayInOut) && (inOutEdge != null))
        {
            var prevMatrixHandles = UnityEditor.Handles.matrix;
            var prevMatrixGizmos = Gizmos.matrix;

            UnityEditor.Handles.matrix = navMeshDisplay.transform.localToWorldMatrix;
            Gizmos.matrix = navMeshDisplay.transform.localToWorldMatrix;

            foreach (var edge in inOutEdge)
            {
                UnityEditor.Handles.DrawBezier(edge.p0, edge.p1, edge.p0, edge.p1, Color.red, null, 15.0f);

                if (displayInOutNormals)
                {
                    var center = (edge.p0 + edge.p1) * 0.5f;
                    Gizmos.color = Color.cyan;

                    float length = 0.1f;

                    Gizmos.DrawLine(center, center + edge.normal * length);

                    Gizmos.color = Color.magenta;
                    Gizmos.DrawLine(center, center + edge.perp * length);
                }
            }

            UnityEditor.Handles.matrix = prevMatrixHandles;
            Gizmos.matrix = prevMatrixGizmos;
        }


        if ((detectInOutProbeDebugView) && (detectInOutProbeDebugViewEdgeIndex > 0))
        {
            var prevMatrixHandles = UnityEditor.Handles.matrix;
            var prevMatrixGizmos = Gizmos.matrix;

            UnityEditor.Handles.matrix = navMeshDisplay.transform.localToWorldMatrix;
            Gizmos.matrix = navMeshDisplay.transform.localToWorldMatrix;

            Boundary boundary = null;
            Mesh mesh = (sourceMesh) ? (sourceMesh.sharedMesh) : (null);

            if ((boundaryDisplay) && (boundaryDisplay.boundary != null))
                boundary = boundaryDisplay.boundary;
            else if ((simplifiedBoundaryDisplay) && (simplifiedBoundaryDisplay.boundary != null))
                boundary = simplifiedBoundaryDisplay.boundary;
            else
            {
                var topology = new Topology(navMeshDisplay.sharedMesh);
                topology.ComputeTriangleNormals();
                boundary = topology.GetBoundary();
            }

            if (boundary != null)
            {
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

                        if (detectInOutProbeDebugViewEdgeIndex > 0)
                        {
                            if (edgeId != detectInOutProbeDebugViewEdgeIndex) continue;
                        }

                        // Get edge
                        var p1 = polyline[edgeCount];
                        var p2 = polyline[(edgeCount + 1) % polyline.Count];

                        UnityEditor.Handles.DrawBezier(p1, p2, p1, p2, Color.yellow, null, 5.0f);
                    }
                }
            }

            UnityEditor.Handles.matrix = prevMatrixHandles;
            Gizmos.matrix = prevMatrixGizmos;
        }

        if ((debugSamplePoints) && (samplePoints != null))
        {
            var prevMatrixHandles = UnityEditor.Handles.matrix;
            var prevMatrixGizmos = Gizmos.matrix;

            UnityEditor.Handles.matrix = navMeshDisplay.transform.localToWorldMatrix;
            Gizmos.matrix = navMeshDisplay.transform.localToWorldMatrix;

            Gizmos.color = new Color(0.0f, 1.0f, 0.0f, 0.5f);
            foreach (var sample in samplePoints)
            {
                Gizmos.DrawSphere(sample, inOutSampleLength * 0.5f);
            }

            UnityEditor.Handles.matrix = prevMatrixHandles;
            Gizmos.matrix = prevMatrixGizmos;
        }

        if ((debugEdges) && (debugEdgeList != null))
        {
            var prevMatrixHandles = UnityEditor.Handles.matrix;
            var prevMatrixGizmos = Gizmos.matrix;

            UnityEditor.Handles.matrix = navMeshDisplay.transform.localToWorldMatrix;
            Gizmos.matrix = navMeshDisplay.transform.localToWorldMatrix;
            
            foreach (var edge in debugEdgeList)
            {
                if (((debugEdgeId == -1) || (debugEdgeId == edge.edgeId)) &&
                    ((debugSubEdgeId == -1) || (debugSubEdgeId == edge.subEdgeId) || (edge.subEdgeId == -1)))
                {
                    UnityEditor.Handles.DrawBezier(edge.start, edge.end, edge.start, edge.end, edge.color, null, 5.0f);
                }
            }

            UnityEditor.Handles.matrix = prevMatrixHandles;
            Gizmos.matrix = prevMatrixGizmos;
        }
    }
}
