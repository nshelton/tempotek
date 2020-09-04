using System;
using System.Collections;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.Linq;
using UnityEngine;

public enum LineType
{
    correlationWeights,
    weights,
    doubleCorr,
    levels,
    histogram,
    lfo,
    beat,
    bpm,
    rollingAvg,
    fft
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
            case LineType.doubleCorr:
                m_targetArray = m_beatDetector.m_doubleCorrelation;
                break;
            case LineType.correlationWeights:
                m_targetArray = m_beatDetector.m_currentWeights;
                break;
            case LineType.weights:
                m_targetArray = m_beatDetector.m_weights;
                break;
            case LineType.levels:
                m_targetArray = m_beatDetector.m_currentLevels;
                break;
            case LineType.histogram:
                m_targetArray = m_beatDetector.m_bpmHistogram;
                break;
            case LineType.rollingAvg:
                m_targetArray = m_beatDetector.m_rollingCorrAverage;
                break;
            case LineType.fft:
                m_targetArray = m_beatDetector.m_fftMag;
                break;
            case LineType.lfo:
                m_targetArray = m_beatDetector.m_beatLFO;
                break;
        }

        if (m_positions == null || m_positions.Length != m_targetArray.Length)
        {
            m_positions = new Vector3[m_targetArray.Length];
            m_lineRenderer.positionCount = m_positions.Length;
        }

        float maxVal = Mathf.Max(m_targetArray.Max(), Mathf.Abs(m_targetArray.Min()));
        
        if ( maxVal > 0 )
            for (int i =0; i < m_targetArray.Length; i++)
            {
                m_positions[i].x = (float)i / m_targetArray.Length;
                m_positions[i].y =  m_targetArray[i] / maxVal;
 
            }

        m_lineRenderer.SetPositions(m_positions);

    }
}
