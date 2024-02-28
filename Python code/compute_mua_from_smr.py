import neo
import numpy
import matplotlib.pyplot as plt


reader = neo.io.Spike2IO(filename='20180808_cell04_D1_DL_sp_activity_BFl_M1.smr',try_signal_grouping=False)
data = reader.read( lazy=False)[0]

#data.segments[0].analogsignals[0].name #to know the channel name
for i in range(0,len(data.segments[0].analogsignals)):
    if (data.segments[0].analogsignals[i].name == 'Vm'): Vm = np.asarray(data.segments[i].analogsignals[0])
    if (data.segments[0].analogsignals[i].name == 'M1'): motor = np.asarray(data.segments[i].analogsignals[0])
    if (data.segments[0].analogsignals[i].name == 'BFl'): somato = np.asarray(data.segments[i].analogsignals[0])
    if (data.segments[0].analogsignals[i].name == 'V1'): visual = np.asarray(data.segments[i].analogsignals[0])

motor = np.asarray(data.segments[0].analogsignals[1])

mua = SpectralMUA_library.computeMUA(motor1939000)
logMUA = np.log(mua)
logMUA = SpectralMUA_library.smooth(logMUA,16) #80 ms smoothing as in Capone Cerebral Cortex 2020
logMUA = logMUA[8:]
parameters, Pfc_UDs_MUA = SpectralMUA_library.compute_UDs_logMUA(np.asarray(np.vstack((logMUA,logMUA))))
#Note that the above line is extracting peaks from the MUA...
#It can be correct, but I still need to validate it
Pfc_UDs = []
for i in range(0,len(Pfc_UDs_MUA)):
    Pfc_UDs.append(np.repeat(Pfc_UDs_MUA[i],16*5))