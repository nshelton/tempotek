using Lasp;
using System.Collections.Generic;
using UnityEngine;
using CircularBuffer;
using System;
using System.Linq;

public class BeatDetector2 : MonoBehaviour
{
    public int RMS_HISTORY_LENGTH = 512;
    public int MAX_SHIFT = 256;

    [SerializeField] Camera m_camera = null;
    [SerializeField] AudioLevelTracker m_level = null;

    CircularBuffer<float> m_levelBuffer;

    public float m_lastProcessingTime;

    public int Divisions = 1;

    public float[] m_autoCorr;
    public float[] m_autCorrHighpass;
    public float[] m_autCorrHighpassSmooth;
    public float[] m_autCorrHighpassSmoothConfidence;

    public float[] m_currentLevels;

    public float m_smoothing = 0.99f;
    public int m_grooveHighpass = 10;
    public int m_grooveLowpass = 1;

    public float m_confidence = 0;

    int MIN_BPM = 60;
    int MAX_BPM = 180;
    int NUM_IMPORTANT_PEAKS = 20;

    public float Sensitivity
    {
        get { return m_smoothing; }
        set { m_smoothing = value; }
    }

    public class Peak
    {
        public float index;
        public float value;
        public float harmonics;
    }

    [NonSerialized] public List<Peak> m_currentPeaks = new List<Peak>();
    [NonSerialized] public List<Peak> m_bestPeaks = new List<Peak>();
    [NonSerialized] public List<Peak> m_grid = new List<Peak>();

    float phaseOffset = 0;

    private float Fract(float i)
    {
        return i - Mathf.Floor(i);
    }

    public float GetBeat(int multiplier)
    {
         return Fract( (Time.time - phaseOffset) * (BPM / multiplier) / 60f); 
    }

    
    public float BPM;

    void Start()
    {
        m_autoCorr = new float[MAX_SHIFT];
        m_autCorrHighpass = new float[MAX_SHIFT];
        m_autCorrHighpassSmooth = new float[MAX_SHIFT];
        m_autCorrHighpassSmoothConfidence = new float[MAX_SHIFT];

        m_levelBuffer = new CircularBuffer<float>(RMS_HISTORY_LENGTH);
        m_currentLevels = new float[RMS_HISTORY_LENGTH];

        for (int i = 0; i < RMS_HISTORY_LENGTH; i++)
        {
            m_levelBuffer.PushBack(0);
        }
        m_bestPeaks.Clear();
        for (int i = 0; i < NUM_IMPORTANT_PEAKS; i++)
        {
            m_bestPeaks.Add(new Peak());
        }

    }

    // Update is called at 50 FPS
    void FixedUpdate()
    {
        if (m_levelBuffer == null)
        {
            Start();
        }

        var startTime = DateTime.Now;

        m_levelBuffer.PushBack(m_level.normalizedLevel);

        var rms = m_levelBuffer.ToArray();
        Array.Copy(rms, m_currentLevels, rms.Length);

        // m_camera.backgroundColor = Color.white * (0.5f + 0.5f * GetBeat(1));

        // Calculate autocorrelation
        ///////////////////////////////////////////////////////
        for (int s = 0; s < MAX_SHIFT; s++)
        {
            float sad = 0;
            for (int i = RMS_HISTORY_LENGTH - MAX_SHIFT - s; i < RMS_HISTORY_LENGTH - s; i++)
            {
                sad += Mathf.Abs(rms[i] - rms[i + s]);
            }
            m_autoCorr[s] = sad;
        }

        // highpass
        ///////////////////////////////////////////////////////
        int windowSize = m_grooveHighpass * 2 + 1;
        float sum = 0;

        for (int i = 0; i < windowSize; i++)
        {
            sum += m_autoCorr[i];
            m_autCorrHighpassSmoothConfidence[i] = 0;
        }

        for (int i = m_grooveHighpass; i < MAX_SHIFT - m_grooveHighpass; i++)
        {
            sum += m_autoCorr[i + m_grooveHighpass];
            sum -= m_autoCorr[i - m_grooveHighpass];

            float rollingAvg = sum / windowSize;
            m_autCorrHighpass[i] = rollingAvg - m_autoCorr[i];
      //      m_confidence += Mathf.Abs(m_rollingCorrAverage[i]) / data.Length;
        }


        float lowpassSum = 0;
        windowSize = m_grooveLowpass * 2 + 1;

        float minCorr = m_autCorrHighpassSmooth.Min();
        float maxCorr = m_autCorrHighpassSmooth.Max();

        float sumSq = 0;

        for(int i  =0; i < MAX_SHIFT; i++)
        {
           float d = m_autCorrHighpassSmooth[i] - m_autCorrHighpass[i];
            sumSq += d * d;
        }
        sumSq /= MAX_SHIFT;
        float weight = m_smoothing * (Mathf.Pow(sumSq * 0.01f, 4f));
       // UnityEngine.Debug.Log("sumSq:" + weight);


        for (int i = 0; i < m_grooveLowpass; i++)
        {
            lowpassSum += m_autCorrHighpass[i];
            m_autCorrHighpassSmoothConfidence[i] = 0;
        }

        for (int i = m_grooveLowpass; i < MAX_SHIFT- m_grooveLowpass; i++)
        {
            lowpassSum += m_autCorrHighpass[i + m_grooveLowpass];
            lowpassSum -= m_autCorrHighpass[i - m_grooveLowpass];
            //m_autCorrHighpassSmooth[i] = Mathf.Lerp(m_autCorrHighpassSmooth[i], lowpassSum /windowSize, weight);


            float newValue = m_autCorrHighpass[i];
            float currentValue = m_autCorrHighpassSmooth[i];

            float variance =  (currentValue - newValue) * (currentValue - newValue);

            m_autCorrHighpassSmoothConfidence[i] = variance;

            m_autCorrHighpassSmooth[i] = Mathf.Lerp(m_autCorrHighpassSmooth[i], lowpassSum / windowSize, m_smoothing);

        }


        /// Get Peaks
        ///////////////////////////////////////////////////////////////
        ///
        m_currentPeaks.Clear();

        for (int i = 2; i < MAX_SHIFT-1; i++)
        {
            if (m_autCorrHighpassSmooth[i] > m_autCorrHighpassSmooth[i - 1] &&
                m_autCorrHighpassSmooth[i] > m_autCorrHighpassSmooth[i + 1])
            {
                float interpolated = getLagrangeMinimum(i, m_autCorrHighpassSmooth);
                m_currentPeaks.Add(new Peak { index = interpolated, value = m_autCorrHighpassSmooth[i] });
            }

        }

        m_currentPeaks.Sort((x, y) => y.value.CompareTo(x.value));

       for(int i = 0; i < NUM_IMPORTANT_PEAKS; i++)
       {
            if ( m_currentPeaks.Count > i)
            {
                m_bestPeaks[i] = (m_currentPeaks[i]);
            }
            else
            {
                m_bestPeaks[i] = new Peak { index = 0, value = 0, harmonics = 1 };
            }
       }


        AnalyzePeaks(m_bestPeaks);

        m_lastProcessingTime = (float)(DateTime.Now - startTime).TotalMilliseconds;
        if (Input.GetKey(KeyCode.Space))
        {
            phaseOffset = Time.time;
        }

    }

    private void AnalyzePeaks(List<Peak> peaks)
    {
        string str = "";
        float bestPeak = peaks[0].index;

        // Test hypothesis: it is beat

        for (int i = 0; i < peaks.Count; i++)
        {
            peaks[i].harmonics = 0;
        }


        for(int p = 0; p < peaks.Count; p++)
        {
            float main = peaks[p].index;

            for (int i = 1; i <= 8; i++)
            for (int ii = 1; ii <= 8; ii++)
            {
                float ratio = (float)i / (float)ii;
                int index = (int)Mathf.Round(main * ratio);
                if (index> 0 && index < m_autCorrHighpassSmooth.Length)
                {
                        peaks[p].harmonics += m_autCorrHighpassSmooth[index];
                }
            }
           peaks[p].harmonics /= 8*8*8;
        }


        float bestGridPeriod = bestPeak;
        float bestScore = 0;

        for (int i = 2; i < 9; i++)
        {
            if (i %2 ==1||i==6)
                continue;
            float candidate = bestPeak * i;
            float score = HowManyOverlap(candidate);
            if (score > bestScore)
            {
                bestGridPeriod = candidate;
                bestScore = score;
            }
        }

        for (int i = 1; i < 9; i++)
        {
            if (i % 2 == 1 || i == 6)
                continue;
            float candidate = bestPeak / (float) i;
            float score = HowManyOverlap(candidate);
            if (score > bestScore)
            {
                bestGridPeriod = candidate;
                bestScore = score;
            }
        }

        Divisions = Mathf.FloorToInt(bestPeak / bestGridPeriod);

        // Find peaks on grid
        float sumX = 0;
        float sumY = 0;

        for (int i = 0; i <= 16; i++)
        {
            float index = i * bestGridPeriod;
            if (index > MAX_SHIFT)
                index = 0;

            for (int ii = 0; ii < m_bestPeaks.Count; ii++)
            {
                float diff = Mathf.Abs(m_bestPeaks[ii].index - index);
                if (diff < 3)
                {
                    sumX += i * m_bestPeaks[ii].harmonics;
                    sumY += m_bestPeaks[ii].index * m_bestPeaks[ii].harmonics;
                    break;
                }
            }

        }

        float bestFitPeriod = sumY / sumX;

        for(int i = 0; i < 4; i++)
        {
            var bpm = GetBPMFromOffset(bestFitPeriod);

            if (bpm > MAX_BPM)
            {
                bestFitPeriod *= 2;
            }
            else if (bpm < MIN_BPM)
            {
                bestFitPeriod /= 2;
            }
        }

        BPM = GetBPMFromOffset(bestFitPeriod);

        // generate fitted grid
        m_grid.Clear();
        for (int i = 0; i <= 16; i++)
        {
            float index = i * bestFitPeriod;
            if (index > MAX_SHIFT)
                index = 0;
            m_grid.Add(new Peak { index = index, harmonics = 1, value = 1 });
        }

    }

    private float HowManyOverlap(float period)
    {
        float overlap = 0;

        for( int i = 0; i < 16; i++)
        {
            float targetIndex = period * (float)i;

            for (int ii = 0; ii < m_bestPeaks.Count; ii++)
            {
                float diff = Mathf.Abs(m_bestPeaks[ii].index - targetIndex);
                if (diff < 1)
                {
                    overlap+= m_bestPeaks[ii].value + m_bestPeaks[ii].harmonics + 1;
                    break;
                }
            }
        }

      /*
        int m = 0;
        for(float i = 0; i < MAX_SHIFT; i+=period)
        {
            int index = Mathf.FloorToInt(i);
            overlap += m_autCorrHighpassSmooth[index];
            m += 1;
            if (m > 1000)
                break;
        }
        */
        return overlap;
    }

    public float getLagrangeMinimum(int offset, float[] m_values)
    {
        // shouldn't really happen
        if (offset <= 0 || offset >= m_values.Length-1)
        {
            return m_values[offset];
        }

        var a = m_values[offset - 1];
        var b = m_values[offset];
        var c = m_values[offset + 1];

        // analytic min/max of the 2nd order lagrange interpolation polynomial 
        float y = (a - c) / (2 * (a - 2 * b + c));

        return offset + y;
    }

    public float GetBPMFromOffset(float offset)
    {
        return 60f * 50f / offset;
    }

}
