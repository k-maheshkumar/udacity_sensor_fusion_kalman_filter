# Write a program to update your mean and variance
# when given the mean and variance of your belief
# and the mean and variance of your measurement.
# This program will update the parameters of your
# belief function.

def update(mean1, var1, mean2, var2):
    new_mean = (var2 * mean1 + var1 * mean2) / (var1 + var2)
    new_var = 1.0 / (1/var1 + 1/var2)
    return [new_mean, new_var]

def predict(mean1, var1, mean2, var2):
    new_mean = mean1 + mean2
    new_var = var1 + var2
    return [new_mean, new_var]

# print update(10.,8.,13., 2.)
# print predict(10., 4., 12., 4.)

measurements = [5., 6., 7., 9., 10.]
motion = [1., 1., 2., 1., 1.]
measurement_sig = 4.
motion_sig = 2.
mu = 0.
sig = 10000.

for ms, mt in zip(measurements, motion):
    mu, sig = update(mu, sig, ms, measurement_sig)
    print ("update: [{}, {}]".format(mu, sig))
    mu, sig = predict(mu, sig, mt, motion_sig)
    print ("predict: [{}, {}]".format(mu, sig))

print [mu, sig]
