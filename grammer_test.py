from scipy.integrate import trapz # 复合梯形积分

wavelength = []
rstrength = []

with open("JSL-82148-1.exp.xy", "r") as f:
    for i in f.readlines():
       wavelength.append(float(i.split()[0]))
       rstrength.append(float(i.split()[1]))

wavelength = list(reversed(wavelength))
rstrength = list(reversed(rstrength))

a = trapz(rstrength, wavelength)
print(a)