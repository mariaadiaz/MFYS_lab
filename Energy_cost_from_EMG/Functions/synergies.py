import numpy as np
import pandas as pd
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF
from scipy.stats import norm
from sklearn.decomposition import NMF
from Functions.emg_preprocess import methods_envelope
from scipy.signal import resample, welch, get_window, find_peaks
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from sklearn.preprocessing import normalize

####################################
#### FUNCTIONS USED AS REFERENCE ###
####################################
# LINKS
# https://stackoverflow.com/questions/22767695/python-non-negative-matrix-factorization-that-handles-both-zeros-and-missing-dat
# https://github.com/hiroyuki-kasai/NMFLibrary/blob/master/demo.py
# https://github.com/romi1502/NMF-matlab/blob/master/python/NMF.py

class Processingcollector():
    def __init__(self):
        self.mean_power = 'No output yet'
        self.mean_ssv = 'No output yet' 
        self.mean_dist = 'No output yet' 
        self.y = 'No output yet'
        self.init = False
        self.results = []

    def organize_synergies(self,W_norm, W2_norm, H2, n):
        # Reorder the synergies of (W2,H2) so each column best matches the
        # corresponding column in W_norm. Matching is done using the maximum
        # cosine similarity (dot product on already-normalized vectors).
        # Already-assigned columns are zeroed out to prevent reuse.

        v_ssv = np.zeros(n)
        loc = list()
        n_W2 = np.zeros(W2_norm.shape)
        n_H2 = np.zeros(H2.shape)
        for ii in range(W_norm.shape[1]):
            for jj in range(W2_norm.shape[1]):
                v_ssv[jj] = np.dot(W_norm[:,ii],W2_norm[:,jj])
            if len(loc) >= 1:
                v_ssv[loc] = 0

            v = np.where(v_ssv == np.max(v_ssv))
            loc.append(v[0][0])

            n_W2[:,ii] = W2_norm[:,v[0][0]]
            n_H2[ii,:] = H2[v[0][0],:]

        return n_W2, n_H2, loc

    def process_emg(self, df, fs, max_r_tibialis, max_l_tibialis, max_r_femoris, max_l_femoris, max_r_gastromed, max_l_gastromed, max_r_gastrolat, max_l_gastrolat):
        # The Delsys stream interleaves EMG with IMU channels, so every 4th
        # column (starting at 0) is the EMG signal for one of the 8 sensors.
        df_emg = df.iloc[:,[0,4,8,12,16,20,24,28]]
        df_emg = df_emg.dropna()

        # General variables
        fs_emg = fs

        # ========================
        # Muscle activity envelope
        # ========================
        # Build the linear envelope of each muscle's raw EMG. methods_envelope
        # handles band-pass, rectification, and low-pass smoothing.

        t_emg = np.linspace(0,len(df_emg)/fs_emg,len(df_emg))
        print("EMG time:", np.round(np.max(t_emg),2))

        #-------------Tibialis anterior envelopes
        var = df_emg.filter(like='RightTibialis').columns
        var_x = df_emg[var.values[0]]
        r_tibialis_env = methods_envelope(var_x,fs_emg,'usual','spikes')

        var = df_emg.filter(like='LeftTibialis').columns
        var_x = df_emg[var.values[0]]
        l_tibialis_env = methods_envelope(var_x,fs_emg,'usual','spikes')

        #-------------Rectus femoris envelopes
        var = df_emg.filter(like='RightRectus').columns
        var_x = df_emg[var.values[0]]
        r_femoris_env = methods_envelope(var_x,fs_emg,'usual','spikes')

        var = df_emg.filter(like='LeftRectus').columns
        var_x = df_emg[var.values[0]]
        l_femoris_env = methods_envelope(var_x,fs_emg,'usual','spikes')

        #-------------Gastrocnemius medialis envelopes
        var = df_emg.filter(like='RightGastroMed').columns
        var_x = df_emg[var.values[0]]
        r_gastromed_env = methods_envelope(var_x,fs_emg,'usual','spikes')

        var = df_emg.filter(like='LeftGastroMed').columns
        var_x = df_emg[var.values[0]]
        l_gastromed_env = methods_envelope(var_x,fs_emg,'usual','spikes')

        #-------------Gastrocnemius lateralis envelopes
        var = df_emg.filter(like='RightGastroLat').columns
        var_x = df_emg[var.values[0]]
        r_gastrolat_env = methods_envelope(var_x,fs_emg,'usual','spikes')

        var = df_emg.filter(like='LeftGastroLat').columns
        var_x = df_emg[var.values[0]]
        l_gastrolat_env = methods_envelope(var_x,fs_emg,'usual','spikes')

        # ======================================
        # Max per muscle for normalization same for all muscles but changing the names

        N = 96 # Moving mean window of 96 points, which translates to ~50 ms (Ubeda et al., 2018).

        if max_r_tibialis <= np.max(np.convolve(r_tibialis_env, np.ones(N)/N, mode='valid')):
            max_r_tibialis = np.max(np.convolve(r_tibialis_env, np.ones(N)/N, mode='valid'))
        if max_l_tibialis <= np.max(np.convolve(l_tibialis_env, np.ones(N)/N, mode='valid')):
            max_l_tibialis = np.max(np.convolve(l_tibialis_env, np.ones(N)/N, mode='valid'))
        if max_r_femoris <= np.max(np.convolve(r_femoris_env, np.ones(N)/N, mode='valid')):
            max_r_femoris = np.max(np.convolve(r_femoris_env, np.ones(N)/N, mode='valid'))
        if max_l_femoris <= np.max(np.convolve(l_femoris_env, np.ones(N)/N, mode='valid')):
            max_l_femoris = np.max(np.convolve(l_femoris_env, np.ones(N)/N, mode='valid'))
        if max_r_gastromed <= np.max(np.convolve(r_gastromed_env, np.ones(N)/N, mode='valid')):
            max_r_gastromed = np.max(np.convolve(r_gastromed_env, np.ones(N)/N, mode='valid'))
        if max_l_gastromed <= np.max(np.convolve(l_gastromed_env, np.ones(N)/N, mode='valid')):
            max_l_gastromed = np.max(np.convolve(l_gastromed_env, np.ones(N)/N, mode='valid'))
        if max_r_gastrolat <= np.max(np.convolve(r_gastrolat_env, np.ones(N)/N, mode='valid')):
            max_r_gastrolat = np.max(np.convolve(r_gastrolat_env, np.ones(N)/N, mode='valid'))
        if max_l_gastrolat <= np.max(np.convolve(l_gastrolat_env, np.ones(N)/N, mode='valid')):
            max_l_gastrolat = np.max(np.convolve(l_gastrolat_env, np.ones(N)/N, mode='valid'))

        print("Max values from muscles:", np.round(max_r_tibialis,5), np.round(max_l_tibialis,5), np.round(max_r_femoris,5), np.round(max_l_femoris,5), np.round(max_r_gastromed,5), np.round(max_l_gastromed,5), np.round(max_r_gastrolat,5), np.round(max_l_gastrolat,5))
        # =====================================================================

        # ---- Tibialis anterior envelope
        r_tibialis_norm = r_tibialis_env / max_r_tibialis
        l_tibialis_norm = l_tibialis_env / max_l_tibialis
        # ---- Rectus femoris envelope
        r_femoris_norm = r_femoris_env / max_r_femoris
        l_femoris_norm = l_femoris_env / max_l_femoris
        # ---- Gastrocnemius medialis envelope
        r_gastromed_norm = r_gastromed_env / max_r_gastromed
        l_gastromed_norm = l_gastromed_env / max_l_gastromed
        # ---- Gastrocnemius lateralis envelope
        r_gastrolat_norm = r_gastrolat_env / max_r_gastrolat
        l_gastrolat_norm = l_gastrolat_env / max_l_gastrolat

        # ========================================================
        # Select only the last 20 seconds of the EMG data
        # ========================================================
        # The last 20 s is assume as the steady-state portion of the window; earlier samples include transients from the SF and consider as adaptation period, thus discarded.
        emg_data_norm = pd.DataFrame(np.transpose(np.array([r_femoris_norm,r_gastrolat_norm,r_gastromed_norm,r_tibialis_norm,l_femoris_norm,l_gastrolat_norm,l_gastromed_norm,l_tibialis_norm])) , columns = list(df_emg.columns))
        x = 20 * fs_emg
        emg_data_norm = emg_data_norm.tail(int(np.round(x)))

        # ========================================================
        # Downsample the signals (if necessary for any processing)
        # ========================================================
        factor = int(np.round(len(t_emg)/3))
        r_tibialis = resample(emg_data_norm['RightTibialisAnterior EMG'].values,factor)
        l_tibialis = resample(emg_data_norm['LeftTibialisAnterior EMG'].values,factor)
        r_femoris = resample(emg_data_norm['RightRectusFemoris EMG'].values,factor)
        l_femoris = resample(emg_data_norm['LeftRectusFemoris EMG'].values,factor)
        r_gastrolat = resample(emg_data_norm['RightGastroLateralis EMG'].values,factor)
        l_gastrolat = resample(emg_data_norm['LeftGastroLateralis EMG'].values,factor)
        r_gastromed = resample(emg_data_norm['RightGastroMedialis EMG'].values,factor)
        l_gastromed = resample(emg_data_norm['LeftGastroMedialis EMG'].values,factor)
        t_emg_new = np.linspace(0,max(t_emg),len(r_tibialis))

        # ======================
        # Calculate muscle power
        # ======================
        # Per-muscle signal power via Welch PSD integrated across frequency.
        # 60% overlap gives a stable estimate over the short 20 s segment.

        win_length = 1; # seconds; freq. resolution 2/Tw=0.4 Hz
        nrsamp = np.floor(win_length*fs_emg)
        overlap = np.floor(0.6*nrsamp)
        the_window = get_window('hann', int(nrsamp))

        emg_preprocessed = emg_data_norm.to_numpy()

        fxx, Pxx = welch(emg_preprocessed, fs=fs_emg, window = the_window, noverlap=overlap, nfft = nrsamp, scaling='density', detrend=False, axis = 0)
        power = np.sum(Pxx, axis = 0) * (fxx[2]-fxx[1])

        # ===========================
        # Calculate muscles synergies
        # ===========================

        emg = np.stack((r_femoris, r_gastrolat, r_gastromed, r_tibialis, l_femoris, l_gastrolat, l_gastromed, l_tibialis), axis=1)
        factor = 3
        fix_val = 1000
        t_fix = np.linspace(0,100,fix_val)  # gait cycle resampled to 0-100%

        x = emg
        t = t_emg_new
        # Detect heel-strikes from gastrocnemius medialis bursts; one peak per
        # stride for each leg. Smoothing is wider here (val=160) than for the
        # normalization envelope to suppress within-burst ripple.
        val = 160
        s_rgastro = np.convolve(x[:,2], np.ones(val)/N, mode='same')
        s_lgastro = np.convolve(x[:,6], np.ones(val)/N, mode='same')

        pks_R,_ = find_peaks(s_rgastro,height = np.mean(s_rgastro)*1.2, distance = int(np.round(fs_emg*0.7/factor)), prominence = np.mean(s_rgastro) + (np.std(s_rgastro)/2))
        pks_L,_ = find_peaks(s_lgastro,height = np.mean(s_lgastro)*1.2, distance = int(np.round(fs_emg*0.7/factor)), prominence = np.mean(s_lgastro) + (np.std(s_lgastro)/2))

        # Slice the signal into stride pairs (peak i → peak i+2 = one full
        # right-leg cycle) and time-warp each to 0-100% so cycles can be
        # averaged together regardless of duration.
        jj = 0
        j = 0
        x_sig = np.zeros([len(t_fix),x.shape[1],int(np.floor(pks_R.shape[0]/2))])
        while jj <= len(pks_R)-3:
            c = jj
            starting_point = np.where(t == t[pks_R[c]])
            starting_point = starting_point[0][0]
            end_point = np.where(t == t[pks_R[c+2]])
            end_point = end_point[0][0]
            t1 = np.linspace(0,100,len(t[starting_point:end_point]))
            xvals = x[starting_point:end_point,:]
            set_interp = interp1d(t1,xvals, kind='linear', axis = 0)
            x_sig[:,:,j] = set_interp(t_fix)
            jj = jj + 2
            j = j + 1

        # Average gait cycle: muscles x time(0-100%). Transposed so NMF rows
        # are muscles and columns are timepoints.
        x_total = np.mean(x_sig,axis=2)
        x_total = np.transpose(x_total)

        # NMF requires non-negative input; interpolation can introduce tiny
        # negative values due to floating-point error.
        neg_vals = np.where(x_total <0)
        if len(neg_vals[0] > 0):
            x_total[neg_vals] = 0

        # 4 synergies is the standard count for lower-limb gait analysis.
        # 'nndsvd' init + fixed random_state make results reproducible across
        # iterations so synergies can be compared between min-2 and min-3.
        np.random.seed(42)
        model = NMF(n_components=4, init='nndsvd', max_iter = 200, random_state= 1)
        W = model.fit_transform(x_total)
        H = model.components_
        W_norm = W / np.linalg.norm(W, axis = 0)

        # Variance Accounted For: how well the 4-synergy reconstruction
        # explains the original muscle activations. Typically > 0.9 for gait.
        SST = np.sum(np.sum((np.transpose(x_total) - np.mean(x_total,axis = 1))**2,axis =1))
        SSE = np.sum(np.sum((x_total - np.dot(W,H))**2,axis = 1))
        VAF = 1 - (SSE/SST)

        return power, W_norm, H, VAF, max_r_tibialis, max_l_tibialis, max_r_femoris, max_l_femoris, max_r_gastromed, max_l_gastromed, max_r_gastrolat, max_l_gastrolat

    def compare_muscle_outputs(self,power1, power2, W_norm, W2_norm, H, H2, n):
        # Compare synergies between two windows (min-2 vs min-3) of the same
        # SF iteration to quantify motor variability — the input the EC model
        # regresses on.

        self.mean_power = np.mean([np.sum(power1),np.sum(power2)])

        [n_W2, n_H2, loc] = self.organize_synergies(W_norm, W2_norm, H2, n)

        # SSV: Synergy similarity vector per matched synergy pair (W_norm columns are unit-norm, so the column-wise dot product is cosine similarity).
        SSV = np.sum(W_norm * n_W2, axis = 0)

        self.mean_ssv = np.arccos(np.mean(SSV)) * (180/np.pi)
        # Euclidean distance between matched activation patterns, averaged
        # across the 4 synergies.
        d1 = np.sqrt(np.sum((H - n_H2)**2,axis=1))
        self.mean_dist = np.mean(d1)

        self.y = self.estimation(self.mean_ssv, self.mean_dist, self.mean_power)

        return self.mean_power, self.mean_ssv, self.mean_dist, self.y, n_H2

    def estimation(self,mean_ssv,mean_dist,mean_power):
        # Linear regression mapping motor-variability features to predicted
        # energy cost. Values taken from a previous publication: Díaz, M. A., et al., (2024). Human-in-the-loop optimization of wearable device parameters using an EMG-based objective function.
        
        beta =[1.193,0.008,2.015,3.396] # Matlab model
        X = np.array([mean_ssv,mean_dist,mean_power])
        # Calibration path (self.init=True) passes vectors of features from
        # multiple iterations at once; the normal path passes one scalar each.
        if self.init == True:
            self.y = beta[0] + X[0,:] * beta[1] + X[1,:] * beta[2] + X[2,:] * beta[3]
            self.init = False
        else:
            self.y = beta[0] + X[0] * beta[1] + X[1] * beta[2] + X[2] * beta[3]
        return self.y

    # ---------------------------------------------------------------------------------
    # ---- Bayesian optimization Functions
    #     
    def bayesian_optimization(self, sf_vals, y_vals):
        # Suggest the next stride-frequency (SF) to test. Search space is
        # 75-126% of preferred SF, matching the protocol limits.
        x_range = np.linspace(75,126, 100)
        sample_x = sf_vals
        sample_y = y_vals

        np.random.seed(42)
        # RBF kernel gives smooth interpolation between observed SFs.
        kernel = 1.0 * RBF(length_scale=5, length_scale_bounds=(1e-1, 1e3))
        gp_model = GaussianProcessRegressor(kernel=kernel, n_restarts_optimizer=100)

        gp_model.fit(sample_x.reshape(-1, 1), sample_y)

        y_pred, y_std = gp_model.predict(x_range.reshape(-1, 1), return_std=True)

        # We are minimizing EC, so "best" is the lowest observed y.
        best_idx = np.argmin(sample_y)
        best_x = sample_x[best_idx]
        best_y = sample_y[best_idx]

        # Expected-improvement acquisition: pick the SF that maximizes the
        # expected drop below the current best.
        ei = self.expected_improvement(x_range, gp_model, best_y)
        new_x = x_range[np.argmax(ei)]

        return y_pred, y_std, ei, new_x, x_range

    def expected_improvement(self, x, gp_model, best_y):
        # Standard EI formula for minimization: positive when the GP mean is
        # below best_y or uncertainty is high enough to plausibly beat it.
        y_pred, y_std = gp_model.predict(x.reshape(-1, 1), return_std=True)
        z = (best_y - y_pred) / y_std
        ei = (best_y - y_pred) * norm.cdf(z) + y_std * norm.pdf(z)
        return ei
    
    # ---------------------------------------------------------------------------------
    # ---- Callback Functions
    #     
    def Partial_results_callback(self):
        return [self.mean_power, self.mean_ssv, self.mean_dist, self.y]

    def results_callback(self):
        if self.ready_to_output == True:
            output_variables = self.compare_muscle_outputs()
            self.results.append(output_variables)
            return self.results


