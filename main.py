import matplotlib.pyplot as plt
import tkinter.filedialog as fd
import numpy as np


file = fd.askopenfilename()
inp = open(file,'r')
newname = file.replace('bin', 'txt')
outp = open(newname,'w')

for line in inp:
    a = line.replace(',', '.')
    outp.write(a)

inp.close()
outp.close()

a = np.loadtxt(newname, usecols=(1))
print(a.mean(),'+-', a.std())
print(a)
plt.plot(a)
plt.show()
transformed = np.fft.fft(a)
plt.plot(transformed)
plt.show()
new_t = transformed[:5000]
new_sig = np.fft.ifft(new_t)
print(new_sig.mean(),'+-', new_sig.std())
plt.plot(new_sig)
plt.show()





