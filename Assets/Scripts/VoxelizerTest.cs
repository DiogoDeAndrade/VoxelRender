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
    [ShowIf("shouldMatchSourceGeometry")]
    public bool         mergeCandidates = true;
    [ShowIf("shouldDetectInOut")]
    public bool         mergeEdges;
    private bool        shouldMergeEdges => mergeEdges && shouldDetectInOut;
    [ShowIf("shouldMergeEdges")]
    public float        mergeAngleTolerance = 10;
    [ShowIf("shouldMergeEdges")]
    public float        mergeDistanceTolerance = 0.25f;

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

    class SamplePoint
    {
        public int     sampleId;
        public Vector3 position;
        public Vector3 normal;
        public int     sourceEdgeId;
        public bool    canStart;
        public bool    canEnd;
    }
    struct SampleList
    {
        public bool                isCW;
        public List<SamplePoint>   points;
    }
    private List<SampleList>   samplePoints;



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
        public  Triangle triangle;
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

    MeshOctree meshOctree;

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
                UnityEditor.EditorUtility.DisplayProgressBar("Building...", "Creating mesh octree...", 0.0f);

                meshOctree = sourceMesh.sharedMesh.GetOctree();

                var t1 = stopwatch.ElapsedMilliseconds;
                Debug.Log($"Mesh octree construction time = {t1 - t0} ms");

                UnityEditor.EditorUtility.DisplayProgressBar("Building...", "Creating nav mesh...", 0.1f);

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

                var t2 = stopwatch.ElapsedMilliseconds;
                Debug.Log("Navmesh generation = " + (t2 - t1));

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
                        long dioc0 = stopwatch.ElapsedMilliseconds;
                        DetectInOutContinuous();
                        long dioc1 = stopwatch.ElapsedMilliseconds;
                        Debug.Log($"Detect In Out Time = {dioc1 - dioc0} ms");
                    }

                    UnityEditor.EditorUtility.DisplayProgressBar("Building...", "Matching source geometry...", 0.76f);

                    if (matchSourceGeometry)
                    {
                        long msg0 = stopwatch.ElapsedMilliseconds;
                        MatchSourceGeometry();
                        long msg1= stopwatch.ElapsedMilliseconds;
                        Debug.Log($"Match source geometry = {msg1 - msg0} ms");
                    }

                    if (shouldMergeEdges)
                    {
                        MergeEdges();
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
                                    triangle = mesh.GetTriangle(submeshIndex, triIndex)
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
                                triangle = mesh.GetTriangle(submeshIndex, triIndex)
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


    bool CanReachInfiniteCircular(MeshOctree meshOctree, Vector3 pt, float range, Vector3 mainDir, Vector3 forwardDir)
    {
        if (meshOctree == null) return true;

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

            Triangle triangle = null;
            float    t = float.MaxValue;
            var      d = probeEnd - pt;
            d /= l;
            bool b = meshOctree.Raycast(pt, d, l, ref triangle, ref t);

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
                    triangle = triangle
                });//*/
            }

            if (b) return false;
        }

        return true;
    }

    bool CanReachInfinitePerpendicular(MeshOctree meshOctree, Vector3 pt0, Vector3 pt1, float range, Vector3 forwardDir)
    {
        if (meshOctree == null) return true;

        Vector3 pt = pt0;
        Vector3 ptInc = (pt1 - pt0) / inOutProbeCount;

        for (int i = 0; i < inOutProbeCount; i++)
        {
            Triangle triangle = null;
            float    t = float.MaxValue;

            bool b = meshOctree.Raycast(pt, forwardDir, range, ref triangle, ref t);

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
                    triangle = triangle
                });//*/
            }

            if (b) return false;

            pt += ptInc;
        }

        return true;
    }

    bool CircularRaycast(MeshOctree mesh, Vector3 pt, float range, Vector3 edgeDir, Vector3 perpDir, bool isFirst, bool isLast, int sampleId)
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

            if (mesh != null)
            {
                Triangle triangle = null;
                float    t = float.MaxValue;
                var d = probeEnd - pt;
                d /= l;
                bool b = mesh.Raycast(pt, d, l, ref triangle, ref t);

                if ((sampleId != -1) && (debugSamplePoints) &&
                    ((debugSampleId < 0) || (debugSampleId == sampleId)))
                {
                    debugRays.Add(new DebugRay()
                    {
                        color = Color.cyan,
                        start = pt,
                        dir = d,
                        dist = l,
                        hit = b,
                        t = t,
                        triangle = triangle
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

        samplePoints = new List<SampleList>();

        for (int boundaryIndex = 0; boundaryIndex < boundary.Count; boundaryIndex++)
        {
            var polyline = boundary.Get(boundaryIndex);

            bool isCW = polyline.isCW();

            var currentSamplePointList = new List<SamplePoint>();

            for (int edgeCount = 0; edgeCount < polyline.Count; edgeCount++)
            {
                edgeId++;

                if (edgeId == detectInOutProbeDebugViewEdgeIndex)
                {
                    Debug.Break();
                }

                UnityEditor.EditorUtility.DisplayProgressBar("Building...", $"Sampling in/out ({edgeId} of {totalEdges})", 0.25f + 0.45f * (float)edgeId / totalEdges);

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
                    sampleId++;

                    bool canStart = CircularRaycast(meshOctree, pt, step, edgeDir, Vector3.Cross(edgeDir, nNorm), true, false, sampleId);
                    bool canEnd = CircularRaycast(meshOctree, pt, step, edgeDir, Vector3.Cross(edgeDir, nNorm), false, true, sampleId);

                    currentSamplePointList.Add(new SamplePoint
                    {
                        sampleId = sampleId,
                        position = pt,
                        normal = n,
                        sourceEdgeId = edgeId,
                        canStart = canStart,
                        canEnd = canEnd
                    });

                    pt = pt + edgeDir * inOutSampleLength;
                    n = n + nInc * inOutSampleLength;
                    nNorm = n.normalized;

                    //if (currentSamplePointList.Count == 10) return;
                }
            }

            if (currentSamplePointList.Count > 0)
            {
                samplePoints.Add(new SampleList() { points = currentSamplePointList, isCW = isCW });
            }
        }

        int sampleCount = sampleId;
        sampleId = 0;

        // Sampling done
        // Run through the lists and create candidate edges
        foreach (var sampleList in samplePoints)
        {
            SamplePoint startPoint = null;
            SamplePoint endPoint = null;

            for (int i = 0; i < sampleList.points.Count; i++)
            {
                sampleId++;

                UnityEditor.EditorUtility.DisplayProgressBar("Building...", "Detecting in/out...", 0.7f + 0.05f * sampleId / sampleCount);

                var currentSample = sampleList.points[i];

                if (startPoint != null)
                {
                    if (endPoint != null)
                    {
                        var edgeDir = (endPoint.position - startPoint.position).normalized;
                        var currentDir = (currentSample.position - endPoint.position).normalized;

                        if (Vector3.Dot(currentDir, edgeDir) < angleTolerance)
                        {
                            // There was a shift in direction on this point, add the valid edge
                            AddInOutEdge(meshOctree, startPoint.position, endPoint.position, startPoint.normal, endPoint.normal, sampleList.isCW);

                            // The endpoint could have been a start point as well, so we continue from there
                            startPoint = endPoint;
                            endPoint = null;
                        }
                    }

                    if (currentSample.canStart)
                    {
                        if (currentSample.canEnd)
                        {
                            // This point can continue the edge (can be start and end)
                            endPoint = currentSample;
                        }
                        else
                        {
                            // It could be a start point, but it can't be an end point, so create the edge that was already built
                            // and start a new one
                            if (endPoint != null)
                            {
                                AddInOutEdge(meshOctree, startPoint.position, endPoint.position, startPoint.normal, endPoint.normal, sampleList.isCW);
                            }

                            startPoint = currentSample;
                            endPoint = null;
                        }
                    }
                    else
                    {
                        if (currentSample.canEnd)
                        {
                            // This point can end the edge, so set it as an endpoint
                            endPoint = currentSample;
                        }
                        // Add the edge (which can include the above point or not)
                        if (endPoint != null)
                        {
                            AddInOutEdge(meshOctree, startPoint.position, endPoint.position, startPoint.normal, endPoint.normal, sampleList.isCW);
                        }
                        // Start a new edge
                        startPoint = endPoint = null;
                    }
                }
                else
                {
                    if (currentSample.canStart)
                    {
                        // Start the edge
                        startPoint = currentSample;
                    }
                }
            }
        }
    }

/*    int debugCount = 0;

    public int addedEdgeTest = 0;
    [ReadOnly] public string rejectionCause;*/

    bool AddInOutEdge(MeshOctree mesh, Vector3 p0, Vector3 p1, Vector3 n1, Vector3 n2, bool isCW)
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
                if (!CanReachInfinitePerpendicular(meshOctree, p0, p1, 4 * agentRadius, perp_dir))
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
                if (!CanReachInfiniteCircular(meshOctree, p0, 4 * agentRadius, main_dir, perp_dir))
                {
                    //rejectionCause = "Interior segment (1)";
                    return false;
                }
                main_dir =-main_dir;
                if (!CanReachInfiniteCircular(meshOctree, p1, 4 * agentRadius, main_dir, perp_dir))
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

        inOutEdge.Add(new Edge() { p0 = p0, p1 = p1, normal = (n1 + n2) * 0.5f, perp = perp_dir.normalized });

        return true;
    }

    void MatchSourceGeometry()
    {
        if (meshOctree == null) return;

        var sourceEdges = inOutEdge;
        inOutEdge = new List<Edge>();

        UnityEditor.EditorUtility.DisplayProgressBar("Building...", $"Generating source mesh topology...", 0.77f);

        var topology = new Topology(sourceMesh.sharedMesh);

        UnityEditor.EditorUtility.DisplayProgressBar("Building...", $"Filtering topology...", 0.78f);
        topology.ComputeTriangleNormals();

        float cosSlopeTolerance = Mathf.Cos(Mathf.Deg2Rad * agentMaxSlope);
        float cosDirectionTolerance = Mathf.Cos(directionAngleTolerance * Mathf.Deg2Rad);
        float cosNormalTolerance = Mathf.Cos(normalAngleTolerance * Mathf.Deg2Rad);
        float toleranceRange = 15.0f * agentRadius;

        int edgeId = 0;

        Debug.Log($"{sourceEdges.Count} source edges, {topology.nEdges} edges in original geometry ");

        foreach (var edge in sourceEdges)
        {
            edgeId++;

            UnityEditor.EditorUtility.DisplayProgressBar("Building...", $"Matching source geometry... ({edgeId} of {sourceEdges.Count})", 0.78f + 0.2f * (float)edgeId / sourceEdges.Count);

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

            List<(int, Vector3, int, Vector3)> candidateEdges = new List<(int, Vector3, int, Vector3)>();

            for (int i = 0; i < topology.nEdges; i++)
            {
                var otherEdgeStruct = topology.GetEdgeStruct(i);
                // Check if this is a valid edge to consider
                if (otherEdgeStruct.triangles.Count > 1)
                {
                    Vector3 baseNormal = otherEdgeStruct.triangles[0].normal;
                    bool    isBoundary = true;
                    for (int k = 1; k < otherEdgeStruct.triangles.Count; k++)
                    {
                        if (Vector3.Dot(baseNormal, otherEdgeStruct.triangles[k].normal) > cosSlopeTolerance)
                        {
                            isBoundary = false;
                            break;
                        }
                    }
                    if (!isBoundary) continue;
                }

                var otherEdge = topology.GetEdgeWithIndices(i);
                var e1 = otherEdge.Item2;
                var e2 = otherEdge.Item4;

                if (Vector3.Distance(e1, edgeStart) > Vector3.Distance(e2, edgeStart))
                {
                    // Switch edges
                    e1 = otherEdge.Item4;
                    e2 = otherEdge.Item2;

                    int swap = otherEdge.Item1;
                    otherEdge.Item1 = otherEdge.Item3;
                    otherEdge.Item3 = swap;
                    otherEdge.Item2 = e1;
                    otherEdge.Item4 = e2;
                }

                if ((Vector3.Distance(e1, edgeCenter) > toleranceRange) ||
                    (Vector3.Distance(e2, edgeCenter) > toleranceRange)) continue;

                // Remove edges that exceed a certain vertical threshould (we never want this, I believe, even if the heuristic then 
                // throws it away
                if (Mathf.Abs(edgeCenter.y - (e1.y + e2.y) * 0.5f) > agentStep * 2.0f) continue;

                var eC = (e1 + e2) * 0.5f;

                // Check direction tolerance
                float dp = Mathf.Abs(Vector3.Dot((e2 - e1).normalized, edgeDir));
                if (dp > cosDirectionTolerance)
                {
                    // Check normal tolerance
                    dp = Mathf.Abs(Vector3.Dot(topology.GetEdgeNormal(i), edge.normal));

                    if (dp > cosNormalTolerance)
                    {
                        // Check if center of candidate is in front of the current edge
                        if (Vector3.Dot((eC - edgeCenter).normalized, edge.perp) > 0)
                        {
                            candidateEdges.Add(otherEdge);

                            //if (candidateEdges.Count > 100) break;
                        }
                    }
                }
            }

            // Merge candidates
            if (mergeCandidates)
            {
                Debug.Log($"Candidate Edges = {candidateEdges.Count}");
                MergeCandidateEdges(candidateEdges);
                Debug.Log($"Candidate Edges = {candidateEdges.Count}");
            }

            int subEdgeId = 0;

            foreach (var otherEdge in candidateEdges)
            {
                subEdgeId++;

                var e1 = otherEdge.Item2;
                var e2 = otherEdge.Item4;

                debugEdgeList.Add(new DebugEdge() { color = Color.green, edgeId = edgeId, subEdgeId = subEdgeId, start = e1, end = e2 });

                var eC = (e1 + e2) * 0.5f;
                var eD = (e2 - e1).normalized;

                // Project both edges to the XZ plane
                var edgeStartXZ = edgeStart.x0z();
                var edgeEndXZ = edgeEnd.x0z();
                var perp = edge.perp.x0z();
                var e1XZ = e1.x0z();
                var e2XZ = e2.x0z();
                Vector3 intersection;
                float   tRay;
                float   tLine;
                float   intersectionScore = 0.0f;

                if (Line.Raycast(edgeStartXZ, perp, toleranceRange, e1XZ, e2XZ, out intersection, out tRay, out tLine))
                {
                    var pt = Line.GetClosestPoint(e1, e2, edgeStart + tRay * edge.perp);
                    if (Vector3.Dot(eD, edgeDir) > 0)
                        e1 = pt;
                    else
                        e2 = pt;

                    intersectionScore += 3;
                }
                else
                {
                    // Compute distance fom ray (how close the raycast was of hitting)
                    if (tLine < 0) tLine = Mathf.Abs(tLine) * (e2XZ - e1XZ).magnitude;
                    else if (tLine > 1) tLine = (tLine - 1) * (e2XZ - e1XZ).magnitude;
                    tLine = 1 - Mathf.Clamp01(tLine / agentRadius);

                    intersectionScore += 3 * tLine;
                }
                if (Line.Raycast(edgeEndXZ, perp, toleranceRange, e1XZ, e2XZ, out intersection, out tRay, out tLine))
                {
                    var pt = Line.GetClosestPoint(e1, e2, edgeEnd + tRay * edge.perp);
                    if (Vector3.Dot(eD, edgeDir) > 0)
                        e2 = pt;
                    else
                        e1 = pt;

                    intersectionScore += 3;
                }
                else
                {
                    // Compute distance fom ray (how close the raycast was of hitting)
                    if (tLine < 0) tLine = Mathf.Abs(tLine) * (e2XZ - e1XZ).magnitude;
                    else if (tLine > 1) tLine = (tLine - 1) * (e2XZ - e1XZ).magnitude;
                    tLine = 1 - Mathf.Clamp01(tLine / agentRadius);

                    intersectionScore += 3 * tLine;
                }

/*                if (debugSubEdgeId != -1)
                {
                    if (subEdgeId == debugEdgeId)
                    {
                        debugRays.Add(new DebugRay()
                        {
                            color = Color.cyan,
                            start = edgeStart,
                            dir = perp,
                            dist = toleranceRange,
                            hit = false,
                            t = float.MaxValue,
                            triangle = null
                        });
                        debugRays.Add(new DebugRay()
                        {
                            color = Color.cyan,
                            start = edgeEnd,
                            dir = perp,
                            dist = toleranceRange,
                            hit = false,
                            t = float.MaxValue,
                            triangle = null
                        });
                    }
                }*/

                // The larger the distance, the larger the score
                float distanceScore = Vector3.Dot(eC.x0z() - edgeStartXZ, edge.perp);
                // The larger the difference in Y, the lower the score
                float yScore = -Mathf.Abs(eC.y - edgeCenter.y);
                // Penalize differences in length, until the maximum size
                float lengthScore = Mathf.Min(0, (Vector3.Distance(e1, e2) - edgeLength));
                // The larger the angular difference, the lower the score (we can use the perp with the direction, because we want it 
                // to be larger with the difference, not with the similarity)
                float angleScore = -Mathf.Abs(Vector3.Dot(perp, (e2XZ - e1XZ).normalized)) * 20.0f;

                float edgeScore = intersectionScore + distanceScore + (yScore * 10) + lengthScore + angleScore;

                //Debug.Log($"Score ({subEdgeId})= {edgeScore} / DeltaY = {yScore} / D = {distanceScore} / DiffLength = {lengthScore} / Intersections = {intersectionScore} / Angle = {angleScore}");

                if (score < edgeScore)
                {
                    candidateStart = e1;
                    candidateEnd = e2;
                    candidateNormal = edge.normal;
                    candidatePerp = edge.perp;
                    score = edgeScore;

                    //Debug.Log($"Best score so far: {score}/{subEdgeId}");
                }

                //debugEdgeList.Add(new DebugEdge() { color = Color.yellow, edgeId = edgeId, subEdgeId = -1, start = e1, end = e2 });
            }

            if (score > -float.MaxValue)
            {
                inOutEdge.Add(new Edge() { p0 = candidateStart, p1 = candidateEnd, normal = candidateNormal, perp = candidatePerp });
            }
        }
    }

    void MergeEdges()
    {
        float   cosTolerance = Mathf.Cos(Mathf.Deg2Rad * mergeAngleTolerance);
        float   dist = mergeDistanceTolerance * agentRadius;
        bool    retry = true;

        while (retry)
        {
            retry = false;

            for (int i = 0; i < inOutEdge.Count; i++)
            {
                var currentEdge = inOutEdge[i];

                for (int j = i + 1; j < inOutEdge.Count; j++)
                {
                    var otherEdge = inOutEdge[j];

                    // Check if the angles are different. Compare perpendicular angles, it's the same as direction in this case
                    // If not, skip this candidate
                    if (Vector3.Dot(currentEdge.perp, otherEdge.perp) < cosTolerance) continue;

                    // Check if endpoints are close
                    if (Vector3.Distance(currentEdge.p0, otherEdge.p0) < dist)
                    {
                        // p0 of both edges "match", build merge edge and remove the other edge
                        currentEdge.p0 = otherEdge.p1;
                        retry = true;
                    }
                    else if (Vector3.Distance(currentEdge.p1, otherEdge.p0) < dist)
                    {
                        // p1 of one edge matches p0 of other edge
                        currentEdge.p1 = otherEdge.p1;
                        retry = true;
                    }
                    else if (Vector3.Distance(currentEdge.p0, otherEdge.p1) < dist)
                    {
                        // p0 of one edge matches p1 of other edge
                        currentEdge.p0 = otherEdge.p0;
                        retry = true;
                    }
                    else if (Vector3.Distance(currentEdge.p1, otherEdge.p1) < dist)
                    {
                        // p1 of both edges match
                        currentEdge.p1 = otherEdge.p0;
                        retry = true;
                    }

                    if (retry)
                    {
                        currentEdge.perp = (currentEdge.perp + otherEdge.perp).normalized;
                        currentEdge.normal = (currentEdge.normal + otherEdge.normal).normalized;

                        inOutEdge[i] = currentEdge;
                        inOutEdge.RemoveAt(j);
                        break;
                    }
                }

                if (retry) break;
            }
        }
    }

    void MergeCandidateEdges(List<(int, Vector3, int, Vector3)> candidates)
    {
        float cosTolerance = Mathf.Cos(Mathf.Deg2Rad * 5);
        bool retry = true;

        while (retry)
        {
            retry = false;

            for (int i = 0; i < candidates.Count; i++)
            {
                var (currentI0, currentP0, currentI1, currentP1) = candidates[i];
                var currentDir = (currentP1 - currentP0);
                var currentLength = currentDir.magnitude;
                currentDir /= currentLength;

                for (int j = i + 1; j < candidates.Count; j++)
                {
                    var (otherI0, otherP0, otherI1, otherP1) = candidates[j];

                    // Check if this edge can be merged (check indices)
                    if ((currentI0 != otherI0) &&
                        (currentI0 != otherI1) &&
                        (currentI1 != otherI0) &&
                        (currentI1 != otherI1)) continue;

                    var otherDir = (otherP1 - otherP0);
                    var otherLength = otherDir.magnitude;
                    otherDir /= otherLength;

                    // Check if the angles are different. 
                    // If not, skip this candidate
                    if (Mathf.Abs(Vector3.Dot(currentDir, otherDir)) < cosTolerance) continue;

                    if (currentI0 == otherI0)
                    {
                        currentI0 = otherI1;
                        currentP0 = otherP1;
                        retry = true;
                    }
                    else if (currentI0 == otherI1)
                    {
                        currentI0 = otherI0;
                        currentP0 = otherP0;
                        retry = true;
                    }
                    else if (currentI1 == otherI0)
                    {
                        currentI1 = otherI1;
                        currentP1 = otherP1;
                        retry = true;
                    }
                    else if (currentI1 == otherI1)
                    {
                        currentI1 = otherI0;
                        currentP1 = otherP1;
                        retry = true;
                    }

                    /*float dist = Mathf.Min(currentLength, otherLength) * 0.05f;

                    // Check if endpoints are close
                    if (Vector3.Distance(currentEdge.Item4, otherEdge.Item2) < dist)
                    {
                        // p0 of both edges "match", build merge edge and remove the other edge
                        currentEdge.Item2 = otherEdge.Item4;
                        currentEdge.Item1 = otherEdge.Item3;
                        retry = true;
                    }
                    else if (Vector3.Distance(currentEdge.Item4, otherEdge.Item2) < dist)
                    {
                        // p1 of one edge matches p0 of other edge
                        currentEdge.Item4 = otherEdge.Item4;
                        currentEdge.Item3 = otherEdge.Item3;
                        retry = true;
                    }
                    else if (Vector3.Distance(currentEdge.Item2, otherEdge.Item4) < dist)
                    {
                        // p0 of one edge matches p1 of other edge
                        currentEdge.Item2 = otherEdge.Item2;
                        currentEdge.Item1 = otherEdge.Item1;
                        retry = true;
                    }
                    else if (Vector3.Distance(currentEdge.Item4, otherEdge.Item4) < dist)
                    {
                        // p1 of both edges match
                        currentEdge.Item4 = otherEdge.Item2;
                        currentEdge.Item3 = otherEdge.Item1;
                        retry = true;
                    }*/

                    if (retry)
                    {
                        candidates[i] = (currentI0, currentP0, currentI1, currentP1);
                        candidates.RemoveAt(j);
                        break;
                    }
                }

                if (retry) break;
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
                    Gizmos.DrawSphere(ray.start + ray.dir * ray.t, 0.01f);

                    if (ray.triangle != null)
                    { 
                        Gizmos.DrawLine(ray.triangle.GetVertex(0), ray.triangle.GetVertex(1));
                        Gizmos.DrawLine(ray.triangle.GetVertex(1), ray.triangle.GetVertex(2));
                        Gizmos.DrawLine(ray.triangle.GetVertex(2), ray.triangle.GetVertex(0));
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

            foreach (var sampleList in samplePoints)
            {
                foreach (var sample in sampleList.points)
                {
                    if ((sample.sampleId != debugSampleId) && (debugSampleId != -1)) continue;

                    if (sample.canStart)
                    {
                        if (sample.canEnd) Gizmos.color = new Color(0.0f, 1.0f, 0.0f, 0.5f);
                        else Gizmos.color = new Color(1.0f, 1.0f, 0.0f, 0.5f);
                    }
                    else
                    {
                        if (sample.canEnd) Gizmos.color = new Color(0.0f, 1.0f, 1.0f, 0.5f);
                        else Gizmos.color = new Color(1.0f, 0.0f, 0.0f, 0.5f);
                    }

                    Gizmos.DrawSphere(sample.position, inOutSampleLength * 0.5f);
                }
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

        /*if (meshOctree != null)
        {
            var prevMatrixHandles = UnityEditor.Handles.matrix;
            var prevMatrixGizmos = Gizmos.matrix;

            UnityEditor.Handles.matrix = navMeshDisplay.transform.localToWorldMatrix;
            Gizmos.matrix = navMeshDisplay.transform.localToWorldMatrix;

            Vector3 pt = new Vector3(1.187499f, 0.1081427f, 0.4017856f);
            Vector3 d = new Vector3(0.241813f, 0.0f, 0.9687138f);
            float l = 0.125f;

            //Gizmos.color = Color.yellow;
            //meshOctree.DrawGizmos(2);

            Gizmos.color = Color.green;
            Triangle triangle = null;
            float    t = 0.0f;
            bool b = meshOctree.RaycastWithGizmos(pt, d, l, ref triangle, ref t);
            if (b)
            {
                Gizmos.color = Color.red;
                triangle.DrawGizmo();
            }

            Gizmos.color = Color.yellow;
            Gizmos.DrawLine(pt, pt + d * l);

            UnityEditor.Handles.matrix = prevMatrixHandles;
            Gizmos.matrix = prevMatrixGizmos;
        }*/
    }

    /*[Button("Build Octree")]
    void BuildOctree()
    {
        meshOctree = sourceMesh.sharedMesh.GetOctree(4);
    }*/
}
