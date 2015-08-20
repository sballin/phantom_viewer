import os
import acquire
import process
import matplotlib.pyplot as plt


def dens_slices():
    shot = 1150528015
    time, dens = acquire.time_dens(shot)
    plt.figure()
    plt.plot(time, dens)
    plt.axvline(.60074, c='r')
    plt.axvline(.62916, c='r')
    plt.axvline(.66772, c='r')
    plt.axvspan(.6, .67676, color='#dddddd', lw=0)
    plt.xlabel('Time (s)')
    plt.ylabel('Line average density')
    plt.xlim([.5, 1])
    plt.show()


def dens_gif():
    # Also did for frames 5384:6000
    shot = 1150528017
    time, dens = acquire.time_dens(shot)
    time_phantom = acquire.gpi_series(shot, 'phantom2', 'time')
    time, dens = process.time_crop((time, dens), time_phantom)

    plt.figure()
    plt.plot(time, dens)
    plt.xlim([time_phantom[0], time_phantom[-1]])
    plt.xlabel('Time (s)')
    plt.ylabel('Line average density')
    plt.axvspan(time_phantom[28846], time_phantom[29807], color='#dddddd', lw=0)
    for i in range(28846, 29807):
        v = plt.axvline(time_phantom[i], color='r')
        plt.savefig('outframe_%05d.png' % i)
        v.remove()
    os.system('convert outframe_*.png -layers optimize out.gif')
    os.system('rm outframe_*.png')


dens_gif()
