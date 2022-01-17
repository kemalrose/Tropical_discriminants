



R = QQ[x, y, t]
f = 10^5 * t^2 + 10^5 * t + 100
g = 1/100 * t^2 + 100 * t + 1
I = ideal(x - f, y - g)
eliminate(I, t)




23.3524   16.3524
