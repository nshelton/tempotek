using Lasp;
using System.Collections.Generic;
using UnityEngine;
using CircularBuffer;
using System;

public class BeatDetector : MonoBehaviour
{
    int RMS_HISTORY_LENGTH = 1024;
    int UPSAMPLE = 4;
    int MAX_SHIFT = 512;

    [SerializeField] Camera m_camera = null;
    [SerializeField] AudioLevelTracker m_level = null;
    [SerializeField] SpectrumAnalyzer m_spectrum = null;

    CircularBuffer<float> m_levelBuffer;
    CircularBuffer<float> m_circularBeatLFOBuffer;

    [SerializeField] float m_correlationLerp = 0.99f;
    [SerializeField] float m_fluxLerp = 0.99f;
    [SerializeField] float m_bpmAttentionSpan = 0.999f;

    public float m_lastProcessingTime;

    public float[] m_weights;
    public float[] m_currentWeights;
    public float[] m_doubleCorrelation;
    public float[] m_rollingCorrAverage;
    public float[] m_fftMag;

    public float[] m_currentLevels;
    public float[] m_bpmHistogram;

    public float[] m_lastSpectrum;
    public float[] m_beatLFO;

    int frameNum = 0;
    public float BPM = 0;
    public float rawFreq = 0;
    public float fitFreq = 0;
    public float fitFreqBPM = 0;
    public int m_grooveHighpass = 10;

    public bool m_useFlux = false;
    public float m_confidence = 0;

    int MIN_BPM = 60;
    int MAX_BPM = 180;

    [NonSerialized] public float[] currentPeaks = new float[32];
    [NonSerialized] public float[] currentPeakWeights = new float[32];
    [NonSerialized] public float[] harmonicPeaks = new float[17];

    float phaseOffset = 0;
    public float harmonicThreshold = 1f;
    public float offsetError = 1f;

    public float GetBeat(int multiplier)
    {
         return Mathf.Sin(Mathf.PI * 2 *  (Time.time - phaseOffset) * (BPM / multiplier) / 60f); 
    }

    private float[] m_hammingLUT;

    void Start()
    {
        m_lastSpectrum = new float[m_spectrum.resolution];

        m_weights = new float[MAX_SHIFT]; 
        m_currentWeights = new float[m_weights.Length];
        m_doubleCorrelation = new float[m_weights.Length];
        m_fftMag = new float[m_weights.Length];
        m_rollingCorrAverage = new float[m_weights.Length];

        m_currentLevels = new float[RMS_HISTORY_LENGTH];
        m_hammingLUT = new float[RMS_HISTORY_LENGTH];
        m_beatLFO = new float[RMS_HISTORY_LENGTH];


        m_bpmHistogram = new float[256];
        m_levelBuffer = new CircularBuffer<float>(RMS_HISTORY_LENGTH);
        m_circularBeatLFOBuffer = new CircularBuffer<float>(RMS_HISTORY_LENGTH);

        for(int i =0;i < RMS_HISTORY_LENGTH; i++)
        {
            m_hammingLUT[i] = 0.54f - 0.46f * Mathf.Cos((2 * Mathf.PI * i) / (RMS_HISTORY_LENGTH - 1));
            m_levelBuffer.PushBack(0);
        }
    }

    private void CalculateAutocorrelation(float[] array, ref float[] acorr, int maxShift)
    {
        for (int s = 0; s < maxShift; s++)
        {
            float sad = 0;
            float numSamples = 0;
            for (int i = 0; i < array.Length; i++)
            {
                if (i + s < array.Length)
                {
                    float window = 1;// m_hammingLUT[i];

                    sad += Mathf.Abs(array[i] - array[i + s]) * window;
                    numSamples += window;
                }
            }

            if (numSamples > 0)
            {
                acorr[s] = sad / numSamples;
            }
        }
    }

    float lastFlux;


    private float GetFlux()
    {
        float flux = 0;

        for(int i =0; i < m_spectrum.resolution; i++)
        {
            float powerCurve = 1f - (float)i / (float)m_spectrum.resolution;

             powerCurve = Mathf.Pow(powerCurve, 4f);

            flux += 10f * powerCurve * Mathf.Abs(m_spectrum.spectrumArray[i] - m_lastSpectrum[i]);
            m_lastSpectrum[i] = m_spectrum.spectrumArray[i];
        }
        flux /= m_spectrum.resolution;


        if (flux == 0)
            return lastFlux;

        if (float.IsNaN(flux) || float.IsInfinity(flux))
        {
            flux = 0;
        }
        
        return flux;
    }


    float lerpedFlux;



    // Update is called at 50 FPS
    void FixedUpdate()
    {
        var startTime = DateTime.Now;
        bool useLevel = false;

        if (m_useFlux)
        {
            lerpedFlux = (1f - m_fluxLerp) * GetFlux() + m_fluxLerp * lerpedFlux;
            m_levelBuffer.PushBack(lerpedFlux);
        }
        else
        {
            m_levelBuffer.PushBack(m_level.normalizedLevel);
        }

        var arr = m_levelBuffer.ToArray();
        Array.Copy(arr, m_currentLevels, arr.Length);

        m_circularBeatLFOBuffer.PushBack(GetBeat(1) * m_confidence);

        m_camera.backgroundColor = Color.white * (0.5f + 0.5f * GetBeat(1));
        var beatArr = m_circularBeatLFOBuffer.ToArray();

        Array.Copy(beatArr, m_beatLFO, beatArr.Length);

        CalculateAutocorrelation(arr, ref m_currentWeights, MAX_SHIFT);

        //CalculateAutocorrelation(m_currentWeights, ref m_doubleCorrelation, MAX_SHIFT);

        float alpha = m_correlationLerp;

        for (int  i = 0; i < m_weights.Length; i++)
        {
            m_weights[i] = alpha * m_weights[i] + (1.0f - alpha) * m_currentWeights[i];
        }

        for (int i = 0; i <currentPeaks.Length; i++)
        {
            currentPeaks[i] = 0;
            currentPeakWeights[i] = 0;
        }

        GetPeaks(m_weights);

        GetHarmonicFromPeaks(ref currentPeaks, ref harmonicPeaks);

        GetMaxBPM();

        GetFFT(m_rollingCorrAverage);


        float thisFreq = BPM / 60f;

        if ( Input.GetKeyDown(KeyCode.Space))
        {
            phaseOffset = Time.time;
        }

        if (Input.GetKeyDown(KeyCode.S))
        {
            var ts = DateTime.Now.Millisecond;
            string filename = @"D:\tempotek\" + ts.ToString();
            UnityEngine.Debug.Log("writing " + filename);
            string[] lines = new string[arr.Length];
            for(int i = 0; i< arr.Length; i++)
            {
                lines[i] = arr[i].ToString();
            }
            System.IO.File.WriteAllLines(filename, lines);

        }

        frameNum++;
        m_lastProcessingTime = (float)(DateTime.Now - startTime).TotalMilliseconds;
    }

    float[] swap = new float[512 * 4];

    private void Upsample(float[] input, ref float[] output)
    {
        for (int i = 0; i < input.Length-1; i++)
        {
            float a = input[i];
            float b = input[i + 1];
            for (int j = 0; j < UPSAMPLE; j++)
            {
                float alpha = (float)j / (float)UPSAMPLE;
                swap[i * UPSAMPLE + j] = alpha * input[i+1] + (1f - alpha) * input[i];
            }
        }

        for (int i = 1; i < output.Length - 1; i++)
        {
            output[i] = swap[i - 1] + swap[i] + swap[i + 1];
            output[i] /= 3f;
        }
    }


    float[] m_real;
    float[] m_imag;

    private void GetFFT(float[] data)
    {
        if (m_real == null || m_real.Length != data.Length)
        {
            m_real = new float[data.Length * UPSAMPLE];
            m_imag = new float[data.Length * UPSAMPLE];
        }

        Upsample(data, ref m_real);

        for(int i = 0; i < m_imag.Length; i++)
            m_imag[i] = 0;

        NAudio.Dsp.FastFourierTransform.FFT(true, 8, m_real, m_imag);

        for(int i = 0; i < m_fftMag.Length; i++)
        {
           m_fftMag[i] = Mathf.Sqrt(m_real[i] * m_real[i] + m_imag[i] * m_imag[i]);
          //  m_fftMag[i] = m_real[i];
        }

        int maxIndex = 0;
        float maxValue = 0;
        
        for(int i = 0; i < m_fftMag.Length; i++)
        {
            if (m_fftMag[i] > maxValue)
            {
                maxValue = m_fftMag[i];
                maxIndex = i;
            }
        }

        //UnityEngine.Debug.Log(16f * getLagrangeMinimum(maxIndex, m_fftMag) * (60f * 50f / 1024f ));
        //UnityEngine.Debug.Log(16f * maxIndex * (60f * 50f / 1024f ));

    }

    private int HowManyOverlap(ref float[] peaks, float freq)
    {
        int n = peaks.Length;
        int k = m_weights.Length / (int)freq;
        int overlap = 0;
        for (int i = 0; i < k; i++)
        {
            for (int j = 0; j < n; j++)
            {
                if(Mathf.Abs(peaks[j] - freq * i) < harmonicThreshold)
                {
                    overlap++;
                }
            }
        }

        return overlap;
    }

    public struct PointF
    {
        public float X;
        public float Y;

        public PointF(float x, float y) : this()
        {
            this.X = x;
            this.Y = y;
        }
    }

    public float FindSlope(List<PointF> points )
    {
        //float[] beatWeight = new float[] { 0, 1, 20, 1, 10, 1, 40, 1, 10 };
        float[] beatWeight = new float[] {0, 5, 10, 5, 20, 5, 10, 5, 100 };
        float avg = 0;
        float weight = 0;

        foreach (PointF p in points)
        {
            if ( p.X > 0 && p.Y > 0 && p.X < beatWeight.Length)
            {
                weight += p.X * p.X * beatWeight[(int)p.X];
                avg += (p.Y / p.X) * (p.X * p.X * beatWeight[(int)p.X]);
            }

        }

        return avg / weight;
    }

    private void GetHarmonicFromPeaks(ref float[] peaks, ref float[] harmonics)
    {
        int n = peaks.Length;
        float bestfreq = -1;
        float bestScore = -1;

        for (int i = 0; i < n; i ++)
        {
            if (peaks[i] == 0)
                continue;

            for (int fundamental = 5; fundamental > 0; fundamental--)
            { 
                float freq = peaks[i] / fundamental;
                
                if ( freq > 10)
                {
                    int score = HowManyOverlap(ref peaks, freq);
                    if (score > bestScore)
                    {
                        bestScore = score;
                        bestfreq = freq;
                    }
                }
            }
        }

  
        rawFreq = bestfreq;

        List<PointF> p = new List<PointF>();

        for(int i = 0; i < n; i ++)
        {
            float k = Mathf.Round(peaks[i] / bestfreq);
            if ( Mathf.Abs(k * bestfreq - peaks[i])  < harmonicThreshold)
            {
                p.Add(new PointF(k, peaks[i]));
            }
        }

        float currentFitFreq = FindSlope(p);

        float error = 0;
        float errorNum = 0;
        for (int i = 0; i < p.Count; i++)
        {
            if (p[i].X > 0)
            {
                error += Mathf.Abs(p[i].Y - p[i].X * currentFitFreq);
                errorNum++;
            }

        }
        error = error / errorNum;

        if ( error > 0)
        {
        
            float weight = 1f - 1f / (1f + Mathf.Exp(-error));
            AddInPeak(fitFreq * 2, weight);

            fitFreq = fitFreq * (1f - weight) + currentFitFreq * weight;
        }


        for (int j = 0; j < harmonics.Length; j++)
        {
            harmonics[j] = fitFreq * j;
        }

        fitFreqBPM= GetBPMFromOffset(fitFreq*2);
    }

    private void GetMaxBPM()
    {
        float maxVal = 0;
        float maxBpm = 0;
        for (int i = 0; i < m_bpmHistogram.Length; i ++)
        {
            if ( m_bpmHistogram[i] > maxVal)
            {
                maxVal = m_bpmHistogram[i];
                maxBpm = getLagrangeMinimum(i, m_bpmHistogram);
            }
        }
        BPM = maxBpm;
    }

    int WINDOW_SIZE = 4;

    private void AddInPeak(float peak, float weight)
    {
        float floatOffset = peak;
        if (floatOffset > 10 && floatOffset < m_bpmHistogram.Length - WINDOW_SIZE)
        {
            float bpm = GetBPMFromOffset(floatOffset);
            if (bpm > MIN_BPM && bpm < MAX_BPM)
            {
                for (int d = -WINDOW_SIZE; d <= WINDOW_SIZE; d++)
                {
                    int d_i = (int)bpm + d;
                    float window = 1 - Mathf.Abs((bpm - d_i) / WINDOW_SIZE);
                    m_bpmHistogram[d_i] += window * weight;
                }
            }
        }

        for (int i = 0; i < m_bpmHistogram.Length; i++)
        {
        m_bpmHistogram[i] *= m_bpmAttentionSpan;
        }
    }

    private float Fract(float i)
    {
        return i - Mathf.Floor(i);
    }


    public void GetPeaks(float[] data)
    {

        m_confidence = 0;

        int radius = m_grooveHighpass;
        int windowSize = radius * 2 + 1;
        float sum = 0;

        for (int i = 0; i < windowSize; i++)
        {
            sum += data[i];
            m_rollingCorrAverage[i] =  sum/(float)(i+1) - data[i];
        }

        for (int i = radius; i < data.Length - radius; i++)
        {
            sum += data[i + radius];
            sum -= data[i - radius];

            float rollingAvg = sum/windowSize;
            m_rollingCorrAverage[i] =  rollingAvg - data[i];
            m_confidence += Mathf.Abs(m_rollingCorrAverage[i]) / data.Length;
        }

        //int peakIndex = 1
        // currentPeaks[0] = 0;
        //for (int i = start + 1; i < data.Length - 1; i++)
        
        var peakIndex = 1;
        currentPeaks[0] = 0;
        bool isPeak = false;
        bool wasIncreasing = false;

        for (int i = 1; i < m_rollingCorrAverage.Length-1; i ++)
        {
            if ( data[i-1] > data[i] && data[i+1] > data[i])
            {
                currentPeaks[peakIndex] = getLagrangeMinimum(i, m_rollingCorrAverage);
                currentPeakWeights[peakIndex] = m_rollingCorrAverage[i];

                peakIndex++;

                if (peakIndex == currentPeaks.Length)
                    break;
            }
        }
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
