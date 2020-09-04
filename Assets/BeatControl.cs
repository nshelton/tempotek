using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class BeatControl : MonoBehaviour
{
    private BeatDetector2 m_beat;

    void Update()
    {
        if (m_beat == null)
        {
            m_beat = GameObject.FindObjectOfType<BeatDetector2>();
        }

        float val = (1f - m_beat.GetBeat(1)) * 5 + 5;

        transform.localScale = new Vector3(val, 0.1f, 0.1f);
    }
}
