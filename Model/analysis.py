import numpy as np
import matplotlib.pyplot as plt

class Process:

    def __init__(self):
        self.prop = {}
        self.f = 0
        self.Ps = {}
        self.PsM = {}
        self.Pb = {}
        self.Pfit = {}

    def powerspectrum(self, data):
        self.prop['fN'] = data.prop['rate']/2
        self.prop['dT'] = 1/data.prop['rate']
        self.prop['Nchannels'] = data.prop['Nchannels']
        self.prop['Tmsr'] = data.prop['Tmsr']

        self.f = np.arange(1/data.prop['Tmsr'], self.prop['fN'], 1/data.prop['Tmsr'])
        indP = np.arange(0, len(self.f), 1)

        if data.prop['Nf'] == 1:
            for j in range(data.prop['Nchannels']):
                fft = (np.fft.fft(data.data[j]))*self.prop['dT']
                P = ((fft*np.conj(fft)).real)/data.prop['Tmsr']
                self.Ps[j] = np.take(P, indP)
        else:
            for i in range(data.prop['Nf']):
                for j in range(data.prop['Nchannels']):
                    print(i,j)
                    if i == 0:
                        fft = (np.fft.fft(data.data[j][:,i]))*self.prop['dT']
                        P = ((fft*np.conj(fft)).real)/data.prop['Tmsr']
                        self.Ps[j] = np.take(P, indP)

                    else:
                        fft = (np.fft.fft(data.data[j][:,i]))*self.prop['dT']
                        P = ((fft*np.conj(fft)).real)/data.prop['Tmsr']
                        self.Ps[j] = np.c_[self.Ps[j], np.take(P, indP)]

    def bin_spectrum(self, Nbin):

        for j in range(self.prop['Nchannels']):
            self.PsM[j] = self.Ps[j].mean(1)

        logf = np.log(self.f)
        hist, bin_edges = np.histogram(logf, Nbin)
        delta = (np.diff(bin_edges)/2)[0]
        bin_centers = bin_edges[:-1]+delta

        ind = {}
        for i in range(Nbin):
            ind[i] = np.where((logf >= (bin_centers[i] - delta)) & (logf < (bin_centers[i] + delta)))
        self.fb = np.empty([Nbin, 1])

        for i in range(Nbin):
            self.fb[i] = np.mean(self.f[ind[i]])

        for j in range(self.prop['Nchannels']):
            for i in range(Nbin):
                if i == 0:
                    self.Pb[j] = np.mean(self.PsM[j][ind[i]])
                else:
                    self.Pb[j] = np.r_[self.Pb[j], np.mean(self.PsM[j][ind[i]])]

        self.remove_nan()
        self.prop['NPb'] = 1


    def remove_nan(self):
        for j in range(4):
                inds = np.where(~np.isnan(self.Pb[j][:]))
                self.Pb[j] = self.Pb[j][inds[0]]
        self.fb = self.fb[inds[0]]


    def plot_pb(self):
        fig,ax = plt.subplots(nrows=2, ncols=2, figsize=(14,14))
        ax[0][0].loglog(self.fb,self.Pb[0][:], label='channel 1')
        ax[0][0].set(xlabel="Frequency (Hz)",
                     ylabel="Power (V2/Hz)",);
        ax[0][0].set_xlim([1E-1, 1E4])
        ax[0][0].set_ylim([1E-13, 1E-7])

        ax[0][1].loglog(self.fb,self.Pb[1], label='channel 2')
        ax[0][1].set(xlabel="Frequency (Hz)",
                     ylabel="Power (V2/Hz)");
        ax[0][1].set_xlim([1E-1, 1E4])
        ax[0][1].set_ylim([1E-13, 1E-7])

        ax[1][0].loglog(self.fb,self.Pb[2][:], label='channel 3')
        ax[1][0].set(xlabel="Frequency (Hz)",
                     ylabel="Power (V2/Hz)");
        ax[1][0].set_xlim([1E-1, 1E4])
        ax[1][0].set_ylim([1E-13, 1E-7])

        ax[1][1].loglog(self.fb,self.Pb[3], label='channel 4')
        ax[1][1].set(xlabel="Frequency (Hz)",
                     ylabel="Power (V2/Hz)");
        ax[1][1].set_xlim([1E-1, 1E4]);
        ax[1][1].set_ylim([1E-13, 1E-7]);

        ax[1][0].loglog(self.fb, self.Pfit[2])
        plt.show()


    def sorensen_fit(self, channel = 2, fit_range = np.array([90, 4000])):

        inds = np.where(((self.fb >= fit_range[0]) & (self.fb < 100)) | ((self.fb <= fit_range[1]) & (self.fb > 400)) | ((self.fb < 360) & (self.fb > 280)))
        fb_fit = self.fb[inds[0]]
        Pb_fit = self.Pb[channel][inds[0]]

        S = {}
        for i in range(3):
            S[i] = {}
        for p in range(3):
            for q in range(3):
                S[p][q] = np.dot(np.transpose(fb_fit**(2*p)), Pb_fit**q)

        a  = ((S[0][1]*S[2][2]-S[1][1]*S[1][2]) / (S[0][2]*S[2][2]-S[1][2]**2))[0];
        b  = ((S[1][1]*S[0][2]-S[0][1]*S[1][2]) / (S[0][2]*S[2][2]-S[1][2]**2))[0];

        self.Pfit[channel] = 1/(a+b*self.fb**2)

        self.prop['fc'] = np.sqrt(a/b).real
        self.prop['D'] = (1/b)*2*(np.pi**2)

        x = min(self.f)/self.prop['fc']
        y = max(self.f)/self.prop['fc']
        s = np.sqrt(np.pi) * (   (2*y)/(1+y**2) - (2*x)/(1+x**2) + 2*np.arctan((y-x)/(1+x*y)) - (4/(y-x)) * (np.arctan((y-x)/(1+x*y)))**2)**(1/2)
        self.prop['sfc'] = s*self.prop['fc']/np.sqrt(np.pi*self.prop['fc']*self.prop['Tmsr'])

        g = np.sqrt(((2*y)/(1+y**2)-(2*x)/(1+x**2) + 2*np.arctan((y-x)/(1+x*y)))/((1+np.pi/2)*(y-x)))
        self.prop['sD'] = self.prop['D']*np.sqrt((1+np.pi/2)/(np.pi*self.prop['fc']*self.prop['Tmsr']))*g*s

        self.plot_pb()
