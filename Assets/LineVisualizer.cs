using System;
using System.Collections;
using System.Collections.Generic;
using System.Collections.Specialized;
using UnityEngine;

public enum LineType
{
    correlationWeights,
    correlationWeightsLaplacian,
    levels,
    histogram,
    lfo,
    beat,
    bpm
}

public class LineVisualizer : MonoBehaviour
{
    public LineRenderer m_lineRenderer;
    public BeatDetector m_beatDetector;

    private Vector3[] m_positions;

    [Range(0,5)]
    public float m_pow = 1.0f;

    private float[] m_targetArray;

    public LineType m_lineType = LineType.correlationWeights;


    void Update()
    {
        switch (m_lineType)
        {
            case LineType.correlationWeightsLaplacian:
                m_targetArray = m_beatDetector.m_currentWeightsLaplacian;
                break;
            case LineType.correlationWeights:
                m_targetArray = m_beatDetector.m_currentWeights;
                break;
            case LineType.levels:
                m_targetArray = m_beatDetector.m_currentLevels;
                break;
            case LineType.histogram:
                m_targetArray = m_beatDetector.m_bpmHistogram;
                break;
        }
       


        if (m_positions == null || m_positions.Length != m_targetArray.Length)
        {
            m_positions = new Vector3[m_targetArray.Length];
            m_lineRenderer.positionCount = m_positions.Length;
        }

        for (int i =0; i < m_targetArray.Length; i++)
        {
            m_positions[i].x = (float)i / m_targetArray.Length;
            m_positions[i].y = Mathf.Pow(m_targetArray[i], m_pow);
 
        }

        m_lineRenderer.SetPositions(m_positions);

    }
}
