using System;
using System.Collections;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.Linq;
using UnityEngine;

public enum LineType2
{
    rms,
    autoCorrelation,
    autoCorrelationHighpass,
    autoCorrelationHighpassSmooth,
    confidence,
    flux,
}

public class LineVisualizer2 : MonoBehaviour
{
    public LineRenderer m_lineRenderer;
    public BeatDetector2 m_beatDetector2;

    private Vector3[] m_positions;

    private float[] m_targetArray;
    public LineType2 m_lineType = LineType2.rms;

    void Update()
    {
        switch (m_lineType)
        {
            case LineType2.rms:
                m_targetArray = m_beatDetector2.m_currentLevels;
                break;

            case LineType2.autoCorrelation:
                m_targetArray = m_beatDetector2.m_autoCorr;
                break;

            case LineType2.autoCorrelationHighpass:
                m_targetArray = m_beatDetector2.m_autCorrHighpass;
                break;

            case LineType2.autoCorrelationHighpassSmooth:
                m_targetArray = m_beatDetector2.m_autCorrHighpassSmooth;
                break;

            case LineType2.confidence:
                m_targetArray = m_beatDetector2.m_autCorrHighpassSmoothConfidence;
                break;

            case LineType2.flux:
                m_targetArray = m_beatDetector2.m_flux;
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
