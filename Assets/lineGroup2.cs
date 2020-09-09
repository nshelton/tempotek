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
    public bool Grid = false;

    // Update is called once per frame
    void Update()
    {
        var peakList = Harmonics ? m_beatDetector.m_grid : m_beatDetector.m_bestPeaks;

        float width = m_beatDetector.MAX_SHIFT;
        if (Grid)
        {
            peakList = m_beatDetector.m_alignedGrid;
            width = m_beatDetector.RMS_HISTORY_LENGTH;
        }

        if (m_bestLine != null && peakList != null && peakList.Count > 0)
        {
            m_bestLine.transform.localPosition = new Vector3(peakList[0].index / width, 0, 0);
        }

        if (m_lines.Count != peakList.Count)
        {
            for (int i = 0; i < m_lines.Count; i++)
            {
                Destroy(m_lines[i]);
            }
            m_lines.Clear();

            for (int i = 0; i < peakList.Count; i++)
            {
                m_lines.Add(Instantiate(m_linePrefab, transform));
            }
        }

        for (int i = 0; i < peakList.Count; i++)
        {
            float pos = (float)peakList[i].index / width;
            m_lines[i].transform.localPosition = new Vector3(pos, 0, 0);
            float height = peakList[i].harmonics;
  

            m_lines[i].transform.localScale = new Vector3(1f, height, 1f);
        }
      
    }
}
