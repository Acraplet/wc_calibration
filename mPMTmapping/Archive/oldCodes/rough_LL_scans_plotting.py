import matplotlib.pyplot as plt
plt.style.use(["science", "notebook", "grid"])
A = [10,15,11,12,9,8]
c = [1109.8,6238,1908,2983,1118,3274]
s = [1274821.2, 649017.9, 1059603.6, 908604.2, 1598726.3,2122302.0]

A_s = [100, 10,15,11,12,9,8, 40, 80]
s = [181238.08, 1274821.2, 649017.9, 1059603.6, 908604.2, 1598726.3,2122302.0, 286419.4, 226810.0]


A_as = [100, 10, 15, 30, 20, 22, 18, 5]
a_s = [20344.5, 10975.4, 10119.8, 15499.9, 12411.3, 13196.9, 11531.0, 223912.3]
plt.plot(A, c, 'x', markersize = 10, label = 'chi2 scan for correcting FileID19 (true alpha = 10cm - absorption)')

plt.xlabel('Test alpha (cm)')

plt.ylabel(f'$\chi^2$')

plt.title('Likelihood scan after correction')
plt.legend()
plt.show()

plt.plot(A_s, s, 'x', markersize = 10, label = 'chi2 scan for correcting FileID16 (true alpha = 10cm - scattering)')

plt.xlabel('Test alpha (cm)')

plt.ylabel(f'$\chi^2$')

plt.title('Likelihood scan after correction')
plt.legend()
plt.show()

plt.plot(A_as, a_s, 'x', markersize = 10, label = 'chi2 scan for correcting FileID16 (true alpha_abs = 20cm, true alpha_scat = 20cm)')

plt.xlabel('Test alpha (cm)')

plt.ylabel(f'$\chi^2$')

plt.title('Likelihood scan after correction')
plt.legend()
plt.show()
