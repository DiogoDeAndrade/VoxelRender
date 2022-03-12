using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using NaughtyAttributes;
using UnityEngine.AI;
using Unity.AI.Navigation;
using static RcdtcsUnityUtils;
using static Recast;

public class TestRecastPackage : MonoBehaviour
{
    [SerializeField] private float agentHeight = 0.2f;
    [SerializeField] private float agentRadius = 0.1f;
    [SerializeField] private float agentClimb = 0.025f;
    [SerializeField] private float cellSize = 0.05f;
    [SerializeField] private float cellHeight = 0.05f;
    [SerializeField] private float regionMinSize = 0.0f;
    [SerializeField] private float regionMergeSize = 5.0f;
    [SerializeField] private bool monotonePartition = false;
    [SerializeField] private float edgeMaxLen = 2.0f;
    [SerializeField] private float edgeMaxError = 5.0f;
    [SerializeField] private int vertsPerPoly = 6;
    [SerializeField] private float detailSampleDist = 6.0f;
    [SerializeField] private float detailMaxError = 1.0f;

    private SystemHelper recast;

    [Button("Build NavMesh")]
    void BuildNavMesh()
    {
        RecastMeshParams navMeshParams = new RecastMeshParams();
        navMeshParams.m_agentHeight = agentHeight;
        navMeshParams.m_agentRadius = agentRadius;
        navMeshParams.m_agentMaxClimb = agentClimb;
        navMeshParams.m_agentMaxSlope = 45;

        navMeshParams.m_cellSize = cellSize;
        navMeshParams.m_cellHeight = cellHeight;

        navMeshParams.m_regionMinSize = regionMinSize;
        navMeshParams.m_regionMergeSize = regionMergeSize;
        navMeshParams.m_monotonePartitioning = monotonePartition;

        navMeshParams.m_edgeMaxLen = edgeMaxLen;
        navMeshParams.m_edgeMaxError = edgeMaxError;
        navMeshParams.m_vertsPerPoly = vertsPerPoly;
        navMeshParams.m_detailSampleDist = detailSampleDist;
        navMeshParams.m_detailSampleMaxError = detailMaxError;

        recast = new SystemHelper();

        recast.SetNavMeshParams(navMeshParams);
        recast.ClearComputedData();
        recast.ClearMesh();

        var meshFilters = GetComponentsInChildren<MeshFilter>();
        foreach (var meshFilter in meshFilters)
        {
            recast.AddMesh(meshFilter.sharedMesh, meshFilter.gameObject);
        }

        recast.ComputeSystem();

        Debug.Log("Build time = " + recast.m_ctx.getAccumulatedTime(Recast.rcTimerLabel.RC_TIMER_TOTAL) + " ms");
    }

    private void OnDrawGizmos()
    {
        if (recast == null) return;

        List<Vector3> poly = new List<Vector3>();

        Vector3 bmin = new Vector3(recast.m_cfg.bmin[0], recast.m_cfg.bmin[1], recast.m_cfg.bmin[2]);

        rcPolyMesh polyMesh = recast.m_pmesh;

        Gizmos.color = Color.yellow;
        for (int i = 0; i < polyMesh.npolys; i++)
        {
            int pIndex = i * polyMesh.nvp * 2;

            poly.Clear();
            for (int j = 0; j < polyMesh.nvp; j++)
            {
                if (polyMesh.polys[pIndex + j] == Recast.RC_MESH_NULL_IDX)
                    break;

                int vIndex = polyMesh.polys[pIndex + j] * 3;

                poly.Add(new Vector3(bmin.x + polyMesh.verts[vIndex + 0] * polyMesh.cs,
                                     bmin.y + polyMesh.verts[vIndex + 1] * polyMesh.ch + 0.001f,
                                     bmin.z + polyMesh.verts[vIndex + 2] * polyMesh.cs));
            }

            for (int j = 0; j < poly.Count; j++)
            {
                Gizmos.DrawLine(poly[j], poly[(j + 1) % poly.Count]);
            }
        }
    }
}
