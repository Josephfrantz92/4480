import sys
import pydub
import pylab
import os
import math
import wave
import struct

def goertzel(samples):
    """
    Implementation of the Goertzel algorithm, useful for calculating individual
    terms of a discrete Fourier transform.

    `samples` is a windowed one-dimensional signal originally sampled at `sample_rate`.

    The function returns 2 arrays, one containing the actual frequencies calculated,
    the second the coefficients `(real part, imag part, power)` for each of those frequencies.
    For simple spectral analysis, the power is usually enough.

    Example of usage :
        
        freqs, results = goertzel(some_samples, 44100, (400, 500), (1000, 1100))
    """
    # We will only be processing audio, so 48kHz is the sample rate we will use.
    # Likewise, since we will only be doing DTMF decoding, we can hard code the
    # frequencies we look for. We take both of these out of the argument list for
    # the algorithm.
    # -Maudrie
    SAMPLE_RATE=48000
    freqs=((697, 770, 852, 941, 1209, 1336, 1477))
    window_size = len(samples)
    f_step = SAMPLE_RATE / float(window_size)
    f_step_normalized = 1.0 / window_size

    # Calculate all the DFT bins we have to compute to include frequencies
    # in `freqs`.
    bins = set()
    for f_range in freqs:
        f_start, f_end = f_range
        k_start = int(math.floor(f_start / f_step))
        k_end = int(math.ceil(f_end / f_step))

        if k_end > window_size - 1: raise ValueError('frequency out of range %s' % k_end)
        bins = bins.union(range(k_start, k_end))

    # For all the bins, calculate the DFT term
    n_range = range(0, window_size)
    freqs = []
    results = []
    for k in bins:

        # Bin frequency and coefficients for the computation
        f = k * f_step_normalized
        w_real = 2.0 * math.cos(2.0 * math.pi * f)
        w_imag = math.sin(2.0 * math.pi * f)

        # Doing the calculation on the whole sample
        d1, d2 = 0.0, 0.0
        for n in n_range:
            y  = samples[n] + w_real * d1 - d2
            d2, d1 = d1, y

        # Storing results `(real part, imag part, power)`
        #We don't need real part or imag part, we only need power
        #-Maudrie
        results.append(
            d2**2 + d1**2 - w_real * d1 * d2
        )
        freqs.append(f * sample_rate)
    return freqs, results




if __name__ == '__main__':

    if (len(sys.argv)>=1):
        inputaudio=sys.argv[1]
    else :
        print("Specify file path as command line arg")

    with contextlib.closing(wave.open(inputaudio,'r')) as f:
        frames = f.getnframes()
        rate = f.getframerate()
        duration = frames / float(rate)
        chunk_no=ceil(duration/.04)

    outputstring=""

    for x in range chunk_no:

        t1=x*.4
        t2=t1+.4
        if (t2>duration):
            t2=duration

        chunk=AudioSegment.from_wav(inputaudio)
        chunk=chunk[t1:t2]
        chunk.export('chunk.wav', format="wav")
        meas_freqs, result= goertzel(chunk)

        highest, second=0, 0
        highest_index, second_index=-1,-1
        for y in range len(result):
            if result[y]>second:
                if result[y]>highest:
                    highest=result[y]
                else:
                    second=result[y]
        
        meas_freq1, meas_freq2= meas_freqs[highest_index], meas_freqs[second_index]            

        actual_freq1=find_most_similar(meas_freq1)
        actual_freq2=find_most_similar(meas_freq2)

        if actual_freq1>actual_freq2:
            actual_freq1, actual_freq2= actual_freq2, actual_freq1

        outputstring=outputstring+dtmf_to_digit(actual_freq1, actual_freq2)

        os.remove("chunk.wav")

    print(outputstring)

def dtmf_to_digit(x,y):
    if x==697:
        if y==1209:
            return "1"
        elif y==1336:
            return "2"
        else:
            return "3"
    elif x==770:
        if y==1209:
            return "4"
        elif y==1336:
            return "5"
        else:
            return "6"
    elif x==852:
        if y==1209:
            return "7"
        elif y==1336:
            return "8"
        else:
            return "9"
    else:
        if y==1209:
            return "*"
        elif y==1336:
            return "0"
        else:
            return "#"

def find_most_similar(meas_freq):

    dmtf_freqs=[697, 770, 852, 941, 1209, 1336, 1477]

    error=10000
    most_similar_index=-1
    for x in len(dtmf_freqs):
        new_error=abs((dtmf_freqs[x]-meas_freq)/meas_freq)
        if new_error<error:
            error=new_error
            most_similar_index=x
    return dtmf_freqs[most_similar_index]