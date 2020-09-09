using Lasp;
using System.Collections.Generic;
using UnityEngine;
using CircularBuffer;
using System;
using System.Linq;
using Microsoft.Win32;

public class BeatDetector2 : MonoBehaviour
{
    public int RMS_HISTORY_LENGTH = 512;
    public int MAX_SHIFT = 256;

    [SerializeField] AudioLevelTracker m_level = null;
    [SerializeField] SpectrumAnalyzer m_spectrum = null;

    CircularBuffer<float> m_levelBuffer;
    CircularBuffer<float> m_fluxBuffer;

    public float m_lastProcessingTime;

    public int Divisions = 1;

    public float[] m_flux;
    public float[] m_corrFFT;
    public float[] m_rawFFT; 

    public float[] m_autoCorr;
    public float[] m_autCorrHighpass;
    public float[] m_autCorrHighpassSmooth;
    public float m_autCorrHighpassSmoothIntegral;
    public float[] m_autCorrHighpassSmoothConfidence;

    public float[] m_currentLevels;

    public float m_smoothing = 0.001f;
    public float m_beatSmooth = 0.001f;
    public float m_phaseSmooth = 0.001f;
    
    public int m_grooveHighpass = 10;
    public int m_grooveLowpass = 1;

    public float m_confidence = 0;
    public float m_bestPeriodSmooth = 0;

    public float[] m_lastSpectrum;

    public float m_learningRate = 4f;
    private float[] m_hammingLUT;

    int MIN_BPM = 60;
    int MAX_BPM = 180;
    int NUM_IMPORTANT_PEAKS = 20;

    public float MATCH_THRESHOLD = 1.5f;

    float m_currentPhaseOffset;

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
    [NonSerialized] public List<Peak> m_alignedGrid = new List<Peak>();

    float phaseOffset = 0;

    private float Fract(float i)
    {
        return i - Mathf.Floor(i);
    }

    public float GetBeat(int multiplier)
    {
        if (m_confidence < -0.1)
            return 0f;
         //return Fract( (Time.time - phaseOffset) * (BPM / Mathf.Pow(Divisions, multiplier-1)) / 60f); 
         return Fract( (Time.time - phaseOffset) * (BPM / multiplier) / 60f); 
    }


    
    public float BPM;


    private float GetFlux()
    {
        float flux = 0;


        for (int i = 0; i < m_spectrum.resolution; i++)
        {
            float powerCurve = 1f / ((float)i +1f);

            flux += powerCurve * Mathf.Abs(m_spectrum.spectrumArray[i] - m_lastSpectrum[i]);
        }

        flux /= m_spectrum.resolution;

        if (float.IsNaN(flux) || float.IsInfinity(flux))
        {
            flux = 0;
        }


        for (int i = 0; i < m_spectrum.resolution; i++)
        {
            if (float.IsNaN(m_lastSpectrum[i]))
                m_lastSpectrum[i] = 0f;

            m_lastSpectrum[i] = Mathf.Lerp(m_spectrum.spectrumArray[i], m_lastSpectrum[i], 0.9f);
        }


        return flux;
    }

    InputStream m_stream;

    void Start()
    {
        m_autoCorr = new float[MAX_SHIFT];
        m_autCorrHighpass = new float[MAX_SHIFT];
        m_autCorrHighpassSmooth = new float[MAX_SHIFT];
        m_autCorrHighpassSmoothConfidence = new float[MAX_SHIFT];

        m_corrFFT = new float[MAX_SHIFT];
        m_rawFFT = new float[MAX_SHIFT];

        m_levelBuffer = new CircularBuffer<float>(RMS_HISTORY_LENGTH);
        m_fluxBuffer = new CircularBuffer<float>(RMS_HISTORY_LENGTH);
        m_currentLevels = new float[RMS_HISTORY_LENGTH];
        m_flux = new float[RMS_HISTORY_LENGTH];

        m_lastSpectrum = new float[m_spectrum.resolution];
        m_hammingLUT = new float[RMS_HISTORY_LENGTH];

        for (int i = 0; i < RMS_HISTORY_LENGTH; i++)
        {
            m_hammingLUT[i] = 0.54f - 0.46f * Mathf.Cos((2 * Mathf.PI * i) / (RMS_HISTORY_LENGTH - 1));
            m_fluxBuffer.PushBack(0);
            m_levelBuffer.PushBack(0);
        }
        m_bestPeaks.Clear();
        for (int i = 0; i < NUM_IMPORTANT_PEAKS; i++)
        {
            m_bestPeaks.Add(new Peak());
        }


        m_stream = AudioSystem.GetDefaultInputStream();

    }
    float fluxLerp = 0;

    // Update is called at 50 FPS
    void FixedUpdate()
    {
        if (m_levelBuffer == null)
        {
            Start();
        }

        var startTime = DateTime.Now;

        m_levelBuffer.PushBack(m_level.normalizedLevel);

        fluxLerp = Mathf.Lerp(fluxLerp, GetFlux(), 0.5f);
        m_fluxBuffer.PushBack(fluxLerp);
        Array.Copy(m_fluxBuffer.ToArray(), m_flux, RMS_HISTORY_LENGTH);

        var rms = m_levelBuffer.ToArray();
        Array.Copy(rms, m_currentLevels, RMS_HISTORY_LENGTH);

        // m_camera.backgroundColor = Color.white * (0.5f + 0.5f * GetBeat(1));

        // Calculate autocorrelation
        ///////////////////////////////////////////////////////
        for (int s = 0; s < MAX_SHIFT; s++)
        {
            float sad = 0;
            for (int i = RMS_HISTORY_LENGTH - MAX_SHIFT - s; i < RMS_HISTORY_LENGTH - s; i++)
            {
                 sad += m_hammingLUT[i] * Mathf.Abs( rms[i] -  rms[i + s]);
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

            //if (i > 0)
            //   m_autCorrHighpass[i] = sum / (float)i - m_autoCorr[i];
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
        m_autCorrHighpassSmoothIntegral = 0;

        for (int i  =0; i < MAX_SHIFT; i++)
        {
           float d = m_autCorrHighpassSmooth[i] - m_autCorrHighpass[i];
            sumSq += d * d;
            m_autCorrHighpassSmoothIntegral += m_autCorrHighpassSmooth[i];
        }
        m_autCorrHighpassSmoothIntegral /= MAX_SHIFT;
        sumSq /= MAX_SHIFT;
        float weight = m_smoothing * (Mathf.Pow(sumSq * 0.01f, 4f));
       // UnityEngine.Debug.Log("sumSq:" + weight);

        for (int i = 0; i < m_grooveLowpass; i++)
        {
            lowpassSum += m_autCorrHighpass[i];
        }

        for (int i = m_grooveLowpass; i < MAX_SHIFT- m_grooveLowpass; i++)
        {
            lowpassSum += m_autCorrHighpass[i + m_grooveLowpass];
            lowpassSum -= m_autCorrHighpass[i - m_grooveLowpass];

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

        //DoFFT();
 

       GetPhaseOffset(m_bestPeriodSmooth, m_flux);

       m_lastProcessingTime = (float)(DateTime.Now - startTime).TotalMilliseconds;
       if (Input.GetKey(KeyCode.Space))
       {
            phaseOffset = Time.time;
       }
    }

    private float[] m_real;
    private float[] m_complex;

    private void DoFFT()
    {
        if (m_real == null || m_real.Length != MAX_SHIFT)
        {
            m_real = new float[MAX_SHIFT];
            m_complex = new float[MAX_SHIFT];
        }

        for(int i = 0; i < MAX_SHIFT; i++)
        {
            m_real[i] = m_autCorrHighpassSmooth[i];
            m_complex[i] = 0f;
        }

        int m = 8;

        NAudio.Dsp.FastFourierTransform.FFT(true, m, m_real, m_complex);

        float maxFFT = 0;
        float maxBin = 0;

        for (int i = 0; i < 128; i++)
        {
            m_corrFFT[i] = m_complex[i];
            if (m_corrFFT[i] > maxFFT)
            {
                maxFFT = m_corrFFT[i];
                maxBin = i;
            }
        }
        maxBin = getLagrangeMinimum((int)maxBin, m_corrFFT);
        Debug.Log(maxBin + " " + (maxBin * 60f * 100f / 256f) + "bpm");
    }

    private void GetPhaseOffset(float period, float[] data)
    {
        if (period < 10 || period > 60)
        {
            return;
        }

        m_currentPhaseOffset++;
        if (m_currentPhaseOffset > period)
            m_currentPhaseOffset -= period;

        float bestScore = -1000f;
        float bestOffset = 0;

        int candidates = Mathf.CeilToInt(period);

        for(int i = 0; i <candidates; i++)
        {
            float score = 0;
            for (int j = 0; j < 32; j ++)
            {
                float floatIndex = data.Length - 1 - i - j*period;
                int index = Mathf.RoundToInt(floatIndex);

                if (index > 0 && index < data.Length)
                {
                    score += data[index];
                }
            }

            if ( score > bestScore)
            {
                bestScore = score;
                bestOffset = i;
            }
        }
     //   UnityEngine.Debug.Log(bestOffset + " \t" + m_currentPhaseOffset);
    //    m_currentPhaseOffset = Mathf.Lerp(m_currentPhaseOffset, bestOffset, m_phaseSmooth);
    //     &&
   // (m_alignedGrid.Count == 0 ||
    //RMS_HISTORY_LENGTH - m_alignedGrid[m_alignedGrid.Count - 1].index >= period)
        if (bestOffset == 0)
        {
            m_alignedGrid.Add(new Peak { index = RMS_HISTORY_LENGTH, harmonics = 1 - bestScore, value = 1 });
        }

        if (m_alignedGrid.Count > 0 && m_alignedGrid[0].index <=0)
        {
            m_alignedGrid.RemoveAt(0);
        }

        for(int i =0; i < m_alignedGrid.Count; i++)
        {
            m_alignedGrid[i].index -= 1;
        }
    }

    private float OffsetScore(float[] data, int offset)
    {
        float sad = 0;

        for(int i = 0; i < data.Length - offset; i++)
        {
            sad += Mathf.Abs(data[i] - data[i + offset]);
        }

        return sad / (data.Length - offset);
    }

    float[] m_candidateIntervals = new float[15] {
        1f/8f , 1f/7f, 1f/6f, 1f/5f, 1f/4f, 1f/3f, 1f/2f, 1f, 2f, 3f, 4f, 5f, 6f, 7f, 8f };

    float[] m_candidateIntervalsEven = new float[7] {
        1f/8f , 1f/4f, 1f/2f, 1f, 2f, 6f, 8f };


    private void CalculateHarmonics(List<Peak> peaks)
    {
        for (int i = 0; i < peaks.Count; i++)
        {
            peaks[i].harmonics = 0;
        }


        for (int p = 0; p < peaks.Count; p++)
        {
            float main = peaks[p].index;

            for (int i = 1; i <= 8; i++)
                for (int ii = 1; ii <= 8; ii++)
                {
                    float ratio = (float)i / (float)ii;
                    int index = (int)Mathf.Round(main * ratio);
                    if (index > 0 && index < m_autCorrHighpassSmooth.Length)
                    {
                        peaks[p].harmonics += m_autCorrHighpassSmooth[index];
                    }
                }
            peaks[p].harmonics /= 8 * 8 * 8;
        }
    }


    private float GetGrid(int i)
    {
        if (i > 8)
            i -= 8;

        switch (i)
        {
            case 0: return 2;
            case 1: return 0.5f;
            case 2: return 1;
            case 3: return 0.5f;
            case 4: return 2;
            case 5: return 0.5f;
            case 6: return 1;
            case 7: return 0.5f;
            case 8: return 2f;
            default: return 1;
        }
    }
    private void AnalyzePeaks2(List<Peak> peaks)
    {

        CalculateHarmonics(peaks);

        float bestPeak = peaks[0].index;
        float bestScore = 100f;
        float bestOffset = 0;
        // find lowest score
        for(int i = 0; i < m_candidateIntervalsEven.Length; i++)
        {
            int offs = Mathf.RoundToInt(m_candidateIntervalsEven[i] * bestPeak);
            if ( offs > 10 && offs < 50)
            {
                float score = OffsetScore(m_autCorrHighpassSmooth, offs);
                if ( score < bestScore)
                {
                    bestScore = score;
                    bestOffset = offs;
                }
            }
        }

        if(bestOffset > 25)
        {
            bestOffset /= 2f;
        }

        if (bestOffset < 15)
        {
            bestOffset *= 2f;
        }

        bestOffset = FitIntervalToPeaks(peaks, bestOffset);

        m_bestPeriodSmooth = Mathf.Lerp(m_bestPeriodSmooth, bestOffset, m_beatSmooth);
        BPM = GetBPMFromOffset(m_bestPeriodSmooth);

        m_confidence = 0;
        float maxCorr = m_autCorrHighpassSmooth.Max();
        float total = 0;
        // generate fitted grid
        m_grid.Clear();
        for (int i = 0; i <= 16; i++)
        {
            float index = i * m_bestPeriodSmooth;
            if (Mathf.RoundToInt(index) < MAX_SHIFT)
            {
                m_grid.Add(new Peak { index = index, harmonics = GetGrid(i), value = 1 });
                float val = m_autCorrHighpassSmooth[Mathf.RoundToInt(index)];
                {
                    m_confidence += val / maxCorr;
                    total++;
                }
            }
        }

        m_confidence /= total / 2f;
        if (float.IsNaN(m_confidence))
        {
            m_confidence = -10f;
        }
      //  m_smoothing = 0.01f * Mathf.Exp(m_learningRate * m_confidence);
      //  m_beatSmooth = 0.01f * Mathf.Exp(m_learningRate * m_confidence);
    }

    private float FitIntervalToPeaks(List<Peak> peaks, float hypothesis)
    {
        // Find peaks on grid
        float sumX = 0;
        float sumY = 0;

        for (int i = 0; i <= 16; i++)
        {
            float index = i * hypothesis;
            if (index > MAX_SHIFT)
                break;

            for (int ii = 0; ii < peaks.Count; ii++)
            {
                float diff = Mathf.Abs(peaks[ii].index - index);
                if (diff < MATCH_THRESHOLD)
                {
                    float weight = 1f;// peaks[i].value;
                    sumX += i * weight;
                    sumY += peaks[ii].index * weight;
                    break;
                }
            }
        }

        float bestFitPeriod = sumY / (float.Epsilon + sumX);

        return bestFitPeriod;
    }

    private void AnalyzePeaks(List<Peak> peaks)
    {
        float bestPeak = peaks[0].index;

        CalculateHarmonics(peaks);

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
                if (diff < MATCH_THRESHOLD)
                {
                    sumX += i * m_bestPeaks[ii].value;
                    sumY += m_bestPeaks[ii].index * m_bestPeaks[ii].value;
                    break;
                }
            }

        }

        float bestFitPeriod = sumY / ( float.Epsilon + sumX);

        float thisScore = 0;
        int total = 0;
        for (int i = 0; i <= 16; i++)
        {
            int index = (int)Mathf.Round(i * bestFitPeriod) ;
             
            if (index < MAX_SHIFT)
            {
                total++;
                thisScore += m_autCorrHighpassSmooth[index];
            }
        }
        thisScore /= total;

       // UnityEngine.Debug.Log(thisScore);


        for (int i = 0; i < 4; i++)
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

        m_bestPeriodSmooth = Mathf.Lerp(m_bestPeriodSmooth, bestFitPeriod, m_beatSmooth);
        BPM = GetBPMFromOffset(m_bestPeriodSmooth);

        m_confidence = 0;
        float maxCorr = m_autCorrHighpassSmooth.Max();
         total = 0;
        // generate fitted grid
        m_grid.Clear();
        for (int i = 0; i <= 16; i++)
        {
            float index = i * m_bestPeriodSmooth;
            if (Mathf.RoundToInt(index) < MAX_SHIFT)
            {
                m_grid.Add(new Peak { index = index, harmonics = GetGrid(i), value = 1 });
                float val = m_autCorrHighpassSmooth[Mathf.RoundToInt(index)];
                {
                    m_confidence += val / maxCorr;
                    total++;
                }
            }
        }

        m_confidence /= total / 2f;
        if (float.IsNaN(m_confidence))
        {
            m_confidence = -10f;
        }
        //  m_smoothing = 0.01f * Mathf.Exp(m_learningRate * m_confidence);
        //  m_beatSmooth = 0.01f * Mathf.Exp(m_learningRate * m_confidence);

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
                if (diff < MATCH_THRESHOLD)
                {
                    overlap += m_bestPeaks[ii].value;
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
