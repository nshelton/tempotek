using System;
using System.Collections;
using System.Collections.Generic;
using System.Collections.Specialized;
using System.Linq;
using UnityEngine;

public class RadialLineVisualizer : MonoBehaviour
{
    public LineRenderer m_lineRenderer;
    public BeatDetector2 m_beatDetector2;

    public float m_frequency = 1f;
    public float m_offset = 1f;

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

 /*           case LineType2.bpmHistogram:
                m_targetArray = m_beatDetector2.m_bestBeatPeriodHistogram;
                break;
 */
        }

        m_frequency =   m_beatDetector2.MAX_SHIFT / (4f *m_beatDetector2.m_bestPeriodSmooth);
        if (m_positions == null || m_positions.Length != m_targetArray.Length)
        {
            m_positions = new Vector3[m_targetArray.Length];
            m_lineRenderer.positionCount = m_positions.Length;
        }

        float maxVal = Mathf.Max(m_targetArray.Max(), Mathf.Abs(m_targetArray.Min()));
        
        if ( maxVal > 0 )
            for (int i =0; i < m_targetArray.Length; i++)
            {
                m_positions[i] = RadialTransform((float)i / (float)m_targetArray.Length, m_targetArray[i] / maxVal);
            }

        m_lineRenderer.SetPositions(m_positions);

    }

    private Vector3 RadialTransform(float phi, float r)
    {
        phi *= m_frequency * Mathf.PI * 2f;
        r += m_offset;

        return new Vector3(
            r * Mathf.Cos(phi),
            r * Mathf.Sin(phi), 
            0);
    }
}
