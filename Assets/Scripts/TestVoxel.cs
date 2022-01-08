using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using NaughtyAttributes;

[RequireComponent(typeof(MeshFilter))]
[RequireComponent(typeof(MeshRenderer))]
public class TestVoxel : MonoBehaviour
{
    public VoxelObject.MaterialMode materialMode = VoxelObject.MaterialMode.Multi;
    [ShowIf("isModeMulti")]
    public Material[]   materialList;

    bool isModeMulti => materialMode == VoxelObject.MaterialMode.Multi;

    [Button("Generate")]
    public void Generate()
    {
        VoxelObject vo = GetComponent<VoxelObject>();
        vo.gridSize = new Vector3Int(9, 9, 9);
        vo.offset = new Vector3(-vo.voxelSize.x * vo.gridSize.x * 0.5f, -vo.voxelSize.y * vo.gridSize.y * 0.5f, -vo.voxelSize.z * vo.gridSize.z * 0.5f);
        vo.materialMode = materialMode;
        vo.AllocateData();
        vo.Set(0, 0, 0, 1);
        vo.Set(8, 0, 0, 2);
        vo.Set(0, 8, 0, 3);
        //vo.Set(8, 8, 0, 4);
        vo.Set(0, 0, 8, 1);
        vo.Set(8, 0, 8, 2);
        vo.Set(0, 8, 8, 3);
        //vo.Set(8, 8, 8, 4);
        for (int i = 0; i < 9; i++)
        {
            vo.Set(i, 4, 4, 2);
            vo.Set(4, i, 4, 3);
            vo.Set(4, 4, i, 5);
        }

        var mesh = vo.GetMesh();

        MeshFilter mf = GetComponent<MeshFilter>();
        mf.mesh = mesh;

        if (isModeMulti)
        {
            List<Material> materials = new List<Material>();
            for (int i = 0; i < 256; i++)
            {
                var matId = vo.GetMaterialId(i);
                if (matId != -1) materials.Add(materialList[materials.Count]);
            }

            MeshRenderer mr = GetComponent<MeshRenderer>();
            mr.materials = materials.ToArray();
        }
    }
}
