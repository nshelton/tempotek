using System.Collections;
using System.Collections.Generic;
using UnityEngine;

public class lineGroup : MonoBehaviour
{

    int numLines = 9;
    public GameObject m_linePrefab;
    public BeatDetector m_beatDetector;

    List<GameObject> m_lines = new List<GameObject>();

    void Start()
    {
        for(int i = 0; i < numLines; i++)
        {
            m_lines.Add(Instantiate(m_linePrefab, transform));
        }
    }

    // Update is called once per frame
    void Update()
    {
        for (int i = 0; i < numLines; i++) {

            m_lines[i].transform.localPosition = new Vector3(m_beatDetector.harmonicPeaks[i], 0, 0);
        }
    }
}
