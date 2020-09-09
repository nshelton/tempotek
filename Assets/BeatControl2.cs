using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class BeatControl2 : MonoBehaviour
{
    private BeatDetector2 m_beat;
    private MeshRenderer m_mesh;

    public int beat = 2;

    void Update()
    {
        if (m_beat == null)
        {
            m_beat = GameObject.FindObjectOfType<BeatDetector2>();
        }
        if (m_mesh == null)
        {
            m_mesh = GetComponent<MeshRenderer>();
        }

        float val = (1f - m_beat.GetBeat(beat));

        
        m_mesh.material.SetColor("_EmissionColor", Color.white * val);
    }
}
