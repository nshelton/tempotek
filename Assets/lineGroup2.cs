using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class lineGroup2 : MonoBehaviour
{

    public GameObject m_linePrefab;
    public BeatDetector2 m_beatDetector;

    List<GameObject> m_lines = new List<GameObject>();
    public GameObject m_bestLine;

    public bool Harmonics = false;
 
    // Update is called once per frame
    void Update()
    {
        var peakList = Harmonics? m_beatDetector.m_grid : m_beatDetector.m_bestPeaks;

        if (m_bestLine != null)
        {
            m_bestLine.transform.localPosition = new Vector3(peakList[0].index / (float)m_beatDetector.MAX_SHIFT, 0, 0);
        }

        if (m_lines.Count != peakList.Count)
        {
            m_lines.Clear();
            for (int i = 0; i < peakList.Count; i++)
            {
                m_lines.Add(Instantiate(m_linePrefab, transform));
            }
        }

        for (int i = 0; i < peakList.Count; i++)
        {
            float pos = (float)peakList[i].index / ((float)m_beatDetector.MAX_SHIFT);
            m_lines[i].transform.localPosition = new Vector3(pos, 0, 0);
            float height = peakList[i].harmonics;
  

            m_lines[i].transform.localScale = new Vector3(1f, height, 1f);
        }
      
    }
}
