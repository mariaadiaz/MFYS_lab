import numpy as np
import matplotlib.pyplot as plt
# import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt, medfilt


# ===================================
# DEFINE FUNCTIONS FOR EMG PROCESSING
# ===================================

def butter_highpass_filter(data, high_cut, fs, order=2):
    """
    Apply Butterworth lowpass filter to data
    :param data: data stream
    :param high_cut: high-cut
    :param fs: sampling rate in Hz
    :param order: order of filter
    :return: y
    """
    b, a = butter_highpass(high_cut, fs, order=order)
    y = filtfilt(b, a, data)
    return y


def butter_highpass(high_cut, fs, order=2):
    """
    Define a Butterworth lowpass
    :param high_cut: high cut-off value
    :param fs: sampling frequency
    :param order: order of filter
    :return: b, a
    """
    nyq = 0.5 * fs
    high = high_cut / nyq
    [b, a] = butter(order, high, btype='high')
    return b, a

def butter_bandpass(low_cut, high_cut, fs, order=2):
    """
    Define a Butterworth bandpass
    :param low_cut: low cut-off value
    :param high_cut: high cut-off value
    :param fs: sampling frequency
    :param order: order of filter
    :return: b, a
    """
    nyq = 0.5 * fs
    low = low_cut / nyq
    high = high_cut / nyq
    [b, a] = butter(order, [low, high], btype='band')
    return b, a


def butter_bandpass_filter(data, low_cut, high_cut, fs, order=2):
    """
    Apply Butterworth bandpass filter to data
    :param data: data stream
    :param low_cut: low-cut
    :param high_cut: high-cut
    :param fs: sampling rate in Hz
    :param order: order of filter
    :return: y
    """
    b, a = butter_bandpass(low_cut, high_cut, fs, order=order)
    y = filtfilt(b, a, data)
    return y


def butter_lowpass_filter(data, high_cut, fs, order=2):
    """
    Apply Butterworth lowpass filter to data
    :param data: data stream
    :param high_cut: high-cut
    :param fs: sampling rate in Hz
    :param order: order of filter
    :return: y
    """
    b, a = butter_lowpass(high_cut, fs, order=order)
    y = filtfilt(b, a, data)
    return y


def butter_lowpass(high_cut, fs, order=2):
    """
    Define a Butterworth lowpass
    :param high_cut: high cut-off value
    :param fs: sampling frequency
    :param order: order of filter
    :return: b, a
    """
    nyq = 0.5 * fs
    high = high_cut / nyq
    [b, a] = butter(order, high, btype='low')
    return b, a

def func_gaussian(t, fwhm):
    return np.exp(-(4 * np.log(2) * t ** 2) / fwhm ** 2)

def sigma2fwhm(sigma):
    return sigma * np.sqrt(8 * np.log(2))


# def medfilt (x, k):
#     """Apply a length-k median filter to a 1D array x.
#     Boundaries are extended by repeating endpoints.
#     """
#     assert k % 2 == 1, "Median filter length must be odd."
#     assert x.ndim == 1, "Input must be one-dimensional."
#     k2 = (k - 1) // 2
#     y = np.zeros ((len (x), k), dtype=x.dtype)
#     y[:,k2] = x
#     for i in range (k2):
#         j = k2 - i
#         y[j:,i] = x[:-j]
#         y[:j,i] = x[0]
#         y[:-j,-(i+1)] = x[j:]
#         y[-j:,-(i+1)] = x[-1]
#     return np.median (y, axis=1)

# ============================================
# Delete spikes from the envelope if necessary
# ============================================

def deletespikes(emg):

    oemg = emg
    filt_emg = medfilt(emg,101)
    max_val = np.max(filt_emg)
    outliers = emg > max_val*1.5
    emg[outliers] = np.nan
    [nans,x] = np.isnan(emg), lambda z: z.nonzero()[0]
    emg[nans]= np.interp(x(nans), x(~nans), emg[~nans])

    plt.figure() 
    plt.plot(oemg, color = 'red')   
    plt.plot(emg,color = 'blue')
    plt.show()

    return emg

# ================
# Methods envelope
# ================
def methods_envelope(emg_signal,fs,method,spikes):

    if method == 'literature':
        data = butter_highpass_filter(signal, 20, fs)
        rect_data = np.abs(data)
        envelope = butter_lowpass_filter(rect_data, 10, fs)

    elif method == 'usual':
        signal = emg_signal - np.mean(emg_signal)
        data = butter_bandpass_filter(signal, 10, 450, fs)
        rect_data = np.abs(data)
        envelope = butter_lowpass_filter(rect_data, 6, fs)

        # plt.figure()
        # plt.plot(emg_signal,color = 'blue')
        # plt.plot(envelope, color = 'red')
        # plt.show()

    if spikes == 'nospikes':
        envelope = deletespikes(envelope)
    
    return envelope