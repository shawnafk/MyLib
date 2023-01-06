sour=np.load("T_sour.npy")
chor=np.load("T_chor.npy")

#asseamble complex scalar
comp_s = sour[...,0] + 1j * sour[...,1]
comp_a = chor[...,0] + 1j * chor[...,1]

#projection, is equivlent to subtract the phase difference from j, we here use the phase subtraction method

phase_chor = np.angle(comp_a)
compensated_sour = comp_s*np.exp(-1j*phase_chor)

plot(imag(compensated_sour[700,:]))
plot(real(compensated_sour[700,:]))
#and the ratio is 
ratio=np.real(compensated_sour)/np.imag(compensated_sour)