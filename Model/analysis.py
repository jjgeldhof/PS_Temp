import numpy as np

class Process:

    def __init__(self):
        self.prop = {}
        self.f = 0
        self.Ps = {}
        self.PsM = {}
        self.Pb = {}

    def powerspectrum(self, data):
        self.prop['fN'] = data.prop['rate']/2
        self.prop['dT'] = 1/data.prop['rate']
        self.prop['Nchannels'] = data.prop['Nchannels']

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

        self.Pb = {}
        for i in range(Nbin):
            for j in range(self.prop['Nchannels']):
                if i == 0:
                    self.Pb[j] = np.mean(self.PsM[j][ind[i]])
                else:
                    self.Pb[j] = np.c_[self.Pb[j], np.mean(self.PsM[j][ind[i]])]
        self.prop['NPb'] = 1
