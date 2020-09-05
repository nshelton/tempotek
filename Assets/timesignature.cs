using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class timesignature : MonoBehaviour
{
    private TextMesh m_text;
    private BeatDetector2 m_beat;

    // Update is called once per frame
    void Update()
    {
        if (m_text == null)
        {
            m_text = GetComponent<TextMesh>();
        }
        if (m_beat == null)
        {
            m_beat = GameObject.FindObjectOfType<BeatDetector2>();
        }

        m_text.text = m_beat.Divisions + " divisions " ;
    }
}
