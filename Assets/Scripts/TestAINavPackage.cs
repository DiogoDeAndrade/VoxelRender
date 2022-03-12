using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using NaughtyAttributes;
using UnityEngine.AI;
using Unity.AI.Navigation;

public class TestAINavPackage : MonoBehaviour
{
    [SerializeField] string agentName;

    void Start()
    {
        
    }

    // Update is called once per frame
    void Update()
    {
        
    }

    [Button("Build NavMesh")]
    void BuildNavMesh()
    {
        NavMeshSurface navmesh = GetComponent<NavMeshSurface>();
        if (navmesh == null)
        {
            navmesh = gameObject.AddComponent<NavMeshSurface>();
        }

        Collider collider = GetComponent<Collider>();
        if (collider != null)
        {
            Bounds localBounds = collider.bounds.ConvertToLocal(transform);

            navmesh.collectObjects = CollectObjects.Volume;
            navmesh.center = localBounds.center;
            navmesh.size = localBounds.size;
        }

        int agentId = -1;
        int settingsCount = NavMesh.GetSettingsCount();
        for (int i = 0; i < settingsCount; i++)
        {
            var settings = NavMesh.GetSettingsByIndex(i);
            var id = settings.agentTypeID;
            var an = NavMesh.GetSettingsNameFromID(id);
            if (an == agentName)
            {
                agentId = id;
            }
        }

        if (agentId != -1)
        {
            navmesh.agentTypeID = agentId;
        }
        else
        {
            Debug.LogError($"Agent {agentName} is not defined!");
        }

        navmesh.BuildNavMesh();
    }
}
