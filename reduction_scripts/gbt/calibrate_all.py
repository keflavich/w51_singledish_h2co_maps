sessions = [14,17,20,21,10,11,16,22]
for sess in sessions:
    execfile("calibrate_session{0}.py".format(sess))
