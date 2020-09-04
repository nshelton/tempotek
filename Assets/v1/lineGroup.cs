using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class lineGroup : MonoBehaviour
{

    public GameObject m_linePrefab;
    public BeatDetector m_beatDetector;

    List<GameObject> m_lines = new List<GameObject>();

    public bool Harmonics = false;
 
    // Update is called once per frame
    void Update()
    {
        var floatList = Harmonics ? m_beatDetector.harmonicPeaks : m_beatDetector.currentPeaks; 
        
        if (m_lines.Count != floatList.Length)
        {
            for (int i = 0; i < floatList.Length; i++)
            {
                m_lines.Add(Instantiate(m_linePrefab, transform));
            }
        }


        for (int i = 0; i < floatList.Length; i++) {
            if ( !float.IsNaN(floatList[i]) && !float.IsInfinity(floatList[i]))
            {
                m_lines[i].transform.localPosition = new Vector3(floatList[i], 0, 0);
            }
        }
    }
}
