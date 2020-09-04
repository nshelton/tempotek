using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class debugText : MonoBehaviour
{

    public BeatDetector m_beat;
    public TextMesh m_text;


    // Update is called once per frame
    void Update()
    {
        m_text.text = m_beat.BPM.ToString("F2") + " \tbpm\n"; 
        m_text.text += (30 * m_beat.m_confidence).ToString("F2") + "\t confidence"; 
    }
}
